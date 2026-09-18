# -*- coding: utf-8 -*-
from Bio.PDB.vectors import rotaxis, Vector
from Bio.PDB import Structure, Model, Chain, PDBParser, PDBIO
from gen_segs import gen_segments_multi
from local_math import angle_between_vectors, dihedral_among_vectors, get_normal_vector
from dihedrals_and_axis import dihedrals_and_axis
from pair_DNA_tetrahedron_v2 import TetrahedronCore, find_valid_configuration
from autofill_sequences import process_autofill, check_input_mismatch
import os
import numpy as np
import math
import warnings
from Bio import BiopythonWarning
import time
import sys
import re
warnings.simplefilter('ignore', BiopythonWarning)


# ========== 1. 读取和分配函数 ==========

def read_config_and_sequences(txt_path: str):
    """读取配置和序列（第1行len_hinge，之后4条链，支持多行）"""
    with open(txt_path, 'r') as f:
        raw_lines = [line.strip() for line in f]
    
    if not raw_lines:
        raise ValueError("空文件")
    
    config_line = raw_lines[0]
    len_hinge = 1
    
    if re.match(r'^\d+$', config_line):
        len_hinge = int(config_line)
    elif '=' in config_line or ':' in config_line:
        separator = '=' if '=' in config_line else ':'
        key, val = config_line.split(separator, 1)
        if 'hinge' in key.lower():
            len_hinge = int(val.strip())
        else:
            raw_lines = [''] + raw_lines
    else:
        raw_lines = [''] + raw_lines
    
    sequence_data = raw_lines[1:]
    sequences_blocks = []
    current_block = []
    
    for line in sequence_data:
        if line == '':
            if current_block:
                sequences_blocks.append(''.join(current_block).upper())
                current_block = []
                if len(sequences_blocks) >= 4:
                    break
        else:
            current_block.append(line)
    
    if current_block and len(sequences_blocks) < 4:
        sequences_blocks.append(''.join(current_block).upper())
    
    if len(sequences_blocks) != 4:
        non_empty = [l.upper() for l in sequence_data if l]
        if len(non_empty) == 4:
            sequences_blocks = non_empty
    
    if len(sequences_blocks) != 4:
        raise ValueError(f"需要4条序列，解析到{len(sequences_blocks)}条")
    
    return len_hinge, sequences_blocks[0], sequences_blocks[1], sequences_blocks[2], sequences_blocks[3]

def fix_ter_placement(pdb_path):
    """
    修复 PDB 中 TER 记录位置：确保每条链的 TER 紧跟在该链原子之后，
    而不是被集中到文件末尾（rearrangeDNA_v2 的副作用）。
    """
    chain_atoms = {}
    ter_records = {}
    other_lines = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                cid = line[21]
                chain_atoms.setdefault(cid, []).append(line)
            elif line.startswith('TER'):
                cid = line[21]
                ter_records[cid] = line
            else:
                other_lines.append(line)

    atom_serial = 1
    with open(pdb_path, 'w') as f:
        for cid in sorted(chain_atoms.keys()):
            for line in chain_atoms[cid]:
                new_line = line[:6] + f"{atom_serial:5d}" + line[11:]
                f.write(new_line)
                atom_serial += 1
            if cid in ter_records:
                ter_line = ter_records[cid]
                new_ter = ter_line[:6] + f"{atom_serial:5d}" + ter_line[11:]
                f.write(new_ter)
                atom_serial += 1
        for line in other_lines:
            f.write(line)


def strip_terminal_phosphate(pdb_path):
    """
    直接修改PDB文件：删除每条链第一个残基（residue sequence number = 1）
    中的 P、OP1、OP2 原子（5'端磷酸基团），以避免 gromacs pdb2gmx 报错。
    """
    with open(pdb_path, 'r') as f:
        lines = f.readlines()

    cleaned = []
    removed_count = 0
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            # PDB 固定列格式（0-based index）：
            # 12-15: atom name, 21: chain ID, 22-25: residue sequence number
            atom_name = line[12:16].strip()
            chain_id = line[21:22]
            res_seq_str = line[22:26].strip()
            
            # 只删链 A/B/C/D 中残基序号为 1 的 P / OP1 / OP2
            if chain_id in 'ABCD' and res_seq_str == '1' and atom_name in ('P', 'OP1', 'OP2'):
                removed_count += 1
                continue
        cleaned.append(line)

    with open(pdb_path, 'w') as f:
        f.writelines(cleaned)
    
    if removed_count:
        print(f"[Clean] Removed {removed_count} terminal phosphate atoms from {pdb_path}")

def auto_assign_sequences(seq_a: str, raw_seqs: list, ignore_last_bases: int):
    """自动分配B、C、D（匹配A[0],A[1],A[2]）"""
    L = len(seq_a)
    if L % 3 != 0:
        raise ValueError(f"A链长度{L}不是3的倍数")
    
    seg_len = L // 3
    a_segments_full = [seq_a[i*seg_len:(i+1)*seg_len] for i in range(3)]
    
    if ignore_last_bases > 0:
        a_segments = [s[:-ignore_last_bases] for s in a_segments_full]
    else:
        a_segments = a_segments_full
    
    comp_map = str.maketrans('ATCG', 'TAGC')
    def revc(s): return s.translate(comp_map)[::-1]
    
    candidates = []
    for idx, seq in enumerate(raw_seqs):
        segs_full = [seq[i*seg_len:(i+1)*seg_len] for i in range(3)]
        if ignore_last_bases > 0:
            segs = [s[:-ignore_last_bases] for s in segs_full]
        else:
            segs = segs_full
        candidates.append({
            'idx': idx,
            'full_seq': seq,
            'match_segs': segs
        })
    
    assignment = {}
    used_indices = set()
    
    for a_seg_idx in range(3):
        target_rev = revc(a_segments[a_seg_idx])
        found = None
        
        for cand in candidates:
            if cand['idx'] in used_indices:
                continue
            for other_seg_idx in range(3):
                if cand['match_segs'][other_seg_idx] == target_rev:
                    found = (cand['idx'], other_seg_idx)
                    break
            if found:
                break
        
        if found is None:
            raise ValueError(f"无法为A[{a_seg_idx}]找到互补段")
        
        cand_idx, matched_seg = found
        assignment[a_seg_idx] = {'raw_idx': cand_idx, 'matched_seg_idx': matched_seg}
        used_indices.add(cand_idx)
    
    seqs = {
        'A': seq_a,
        'B': raw_seqs[assignment[0]['raw_idx']],
        'C': raw_seqs[assignment[1]['raw_idx']],
        'D': raw_seqs[assignment[2]['raw_idx']]
    }
    
    return seqs



# ========== 3. 旋转函数 ==========

def Cn_rotate(structure, n: int, axis: Vector, p0: Vector = Vector(0, 0, 0)) -> None:
    """C_n旋转对称操作"""
    if hasattr(structure, '__getitem__') and not isinstance(structure, list):
        structure = list(structure)
    if len(structure) % n:
        raise ValueError('残基数必须能被 n 整除')
    
    total_res = len(structure)
    res_per_seg = total_res // n

    def rot_mat(unit_ax: np.ndarray, theta: float) -> np.ndarray:
        c, s = np.cos(theta), np.sin(theta)
        ic = 1.0 - c
        ux, uy, uz = unit_ax
        return np.array([
            [c + ux*ux*ic,   ux*uy*ic - uz*s, ux*uz*ic + uy*s],
            [uy*ux*ic + uz*s, c + uy*uy*ic,   uy*uz*ic - ux*s],
            [uz*ux*ic - uy*s, uz*uy*ic + ux*s, c + uz*uz*ic]
        ], dtype=float)

    axis_np = np.array(axis.get_array(), dtype=float)
    axis_unit = axis_np / np.linalg.norm(axis_np)
    p0_np = np.array(p0.get_array(), dtype=float)

    for seg in range(n):
        start = seg * res_per_seg
        end = (seg + 1) * res_per_seg if seg != n - 1 else total_res
        R = rot_mat(axis_unit, 2 * np.pi * seg / n)
        for res in structure[start:end]:
            for atom in res.get_atoms():
                coord = np.array(atom.get_coord(), dtype=float)
                shifted = coord - p0_np
                rotated = R @ shifted
                atom.set_coord(rotated + p0_np)


def rotate_head_alpha(residues, theta_deg, n=1):
    theta = -math.radians(theta_deg)
    nres = len(residues)
    if n <= 0 or n > nres - 1:
        raise ValueError("n 必须满足 1 ≤ n ≤ len(residues)-1")
    for k in range(n):
        bend_at = k
        atom_P  = residues[bend_at]['P']
        atom_O5 = residues[bend_at]['O5\'']
        axis_vec = Vector(atom_P.get_coord()) - Vector(atom_O5.get_coord())
        pivot    = Vector(atom_O5.get_coord())
        rot_mat  = rotaxis(theta, axis_vec)
        for atom in residues[bend_at]:
            if atom.get_name() not in ['OP1', 'OP2']:
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
        for idx in range(bend_at):
            for atom in residues[idx]:
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
    return residues


def rotate_head_gamma(residues, theta_deg, n=1):
    theta = -math.radians(theta_deg)
    nres = len(residues)
    if n <= 0 or n > nres:
        raise ValueError("n 必须满足 1 ≤ n ≤ len(residues)")
    for k in range(n):
        bend_at = k
        atom_C5 = residues[bend_at]['C5\'']
        atom_C4 = residues[bend_at]['C4\'']
        axis_vec = Vector(atom_C5.get_coord()) - Vector(atom_C4.get_coord())
        pivot   = Vector(atom_C4.get_coord())
        rot_mat = rotaxis(theta, axis_vec)
        for atom in residues[bend_at]:
            if atom.get_name() not in ['C5\'', 'O5\'', 'P', 'OP1', 'OP2']:
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
        for idx in range(bend_at):
            for atom in residues[idx]:
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
    return residues


def rotate_head_zeta(residues, theta_deg, n=1):
    theta = -math.radians(theta_deg)
    nres = len(residues)
    if n <= 0 or n > nres - 1:
        raise ValueError("n 必须满足 1 ≤ n ≤ len(residues)-1")
    for k in range(n):
        bend_at = k + 1
        atom_O3 = residues[bend_at]['O3\'']
        atom_P  = residues[bend_at - 1]['P']
        axis_vec = Vector(atom_O3.get_coord()) - Vector(atom_P.get_coord())
        pivot   = Vector(atom_P.get_coord())
        rot_mat = rotaxis(theta, axis_vec)
        coord = np.array(atom_O3.get_coord())
        new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
        atom_O3.set_coord(new_c)
        for idx in range(bend_at):
            for atom in residues[idx]:
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
    return residues


def rotate_head_epsilon(residues, theta_deg, n=1):
    theta = -math.radians(theta_deg)
    nres = len(residues)
    if n <= 0 or n > nres:
        raise ValueError("n 必须满足 1 ≤ n ≤ len(residues)")
    for k in range(n):
        bend_at = k
        atom_C3 = residues[bend_at]['C3\'']
        atom_O3 = residues[bend_at]['O3\'']
        axis_vec = Vector(atom_C3.get_coord()) - Vector(atom_O3.get_coord())
        pivot   = Vector(atom_O3.get_coord())
        rot_mat = rotaxis(theta, axis_vec)
        coord = np.array(atom_O3.get_coord())
        new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
        atom_O3.set_coord(new_c)
        for idx in range(bend_at):
            for atom in residues[idx]:
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
    return residues


def rotate_alpha(residues, theta_deg, n=1):
    theta = -math.radians(theta_deg)
    nres = len(residues)
    if n <= 0 or n > nres - 1:
        raise ValueError("n 必须满足 1 ≤ n ≤ len(residues)-1")
    for k in range(1, n + 1):
        bend_at = nres - k - 1
        atom_P  = residues[bend_at]['P']
        atom_O5 = residues[bend_at]['O5\'']
        axis_vec = Vector(atom_O5.get_coord()) - Vector(atom_P.get_coord())
        pivot    = Vector(atom_P.get_coord())
        rot_mat  = rotaxis(theta, axis_vec)
        for atom in residues[bend_at]:
            if atom.get_name() in ['OP1', 'OP2']:
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
        for idx in range(bend_at + 1, nres):
            for atom in residues[idx]:
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
    return residues


def rotate_gamma(residues, theta_deg, n=1):
    theta = -math.radians(theta_deg)
    nres = len(residues)
    if n <= 0 or n > nres:
        raise ValueError("n 必须满足 1 ≤ n ≤ len(residues)")
    for k in range(1, n + 1):
        bend_at = nres - k
        atom_C5 = residues[bend_at]['C5\'']
        atom_C4 = residues[bend_at]['C4\'']
        axis_vec = Vector(atom_C4.get_coord()) - Vector(atom_C5.get_coord())
        pivot   = Vector(atom_C5.get_coord())
        rot_mat = rotaxis(theta, axis_vec)
        for atom in residues[bend_at]:
            if atom.get_name() in ['C5\'', 'O5\'', 'P', 'OP1', 'OP2']:
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
        for idx in range(bend_at + 1, nres):
            for atom in residues[idx]:
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
    return residues


def rotate_zeta(residues, theta_deg, n=1):
    theta = -math.radians(theta_deg)
    nres = len(residues)
    if n <= 0 or n > nres - 1:
        raise ValueError("n 必须满足 1 ≤ n ≤ len(residues)-1")
    for k in range(1, n + 1):
        bend_at = nres - k
        atom_O3 = residues[bend_at]['O3\'']
        atom_P  = residues[bend_at - 1]['P']
        axis_vec = Vector(atom_P.get_coord()) - Vector(atom_O3.get_coord())
        pivot   = Vector(atom_O3.get_coord())
        rot_mat = rotaxis(theta, axis_vec)
        for atom in residues[bend_at]:
            if atom.get_name() != 'O3\'':
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
        for idx in range(bend_at + 1, nres):
            for atom in residues[idx]:
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
    return residues


def rotate_epsilon(residues, theta_deg, n=1):
    theta = -math.radians(theta_deg)
    nres = len(residues)
    if n <= 0 or n > nres:
        raise ValueError("n 必须满足 1 ≤ n ≤ len(residues)")
    for k in range(1, n + 1):
        bend_at = nres - k
        atom_C3 = residues[bend_at]['C3\'']
        atom_O3 = residues[bend_at]['O3\'']
        axis_vec = Vector(atom_O3.get_coord()) - Vector(atom_C3.get_coord())
        pivot   = Vector(atom_C3.get_coord())
        rot_mat = rotaxis(theta, axis_vec)
        for atom in residues[bend_at]:
            if atom.get_name() != 'O3\'':
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
        for idx in range(bend_at + 1, nres):
            for atom in residues[idx]:
                coord = np.array(atom.get_coord())
                new_c = rot_mat @ (coord - pivot.get_array()) + pivot.get_array()
                atom.set_coord(new_c)
    return residues


def rotate_around_axis(structure, angle_deg: float, axis_vec, pivot) -> None:
    if not isinstance(axis_vec, Vector):
        axis_vec = Vector(*axis_vec)
    if not isinstance(pivot, Vector):
        pivot = Vector(*pivot)
    theta = np.deg2rad(angle_deg)
    axis_np = np.array(axis_vec.get_array(), dtype=float)
    axis_unit = axis_np / (np.linalg.norm(axis_np) + 1e-16)
    ux, uy, uz = axis_unit
    c, s = np.cos(theta), np.sin(theta)
    ic = 1.0 - c
    R = np.array([
        [c + ux*ux*ic, ux*uy*ic - uz*s, ux*uz*ic + uy*s],
        [uy*ux*ic + uz*s, c + uy*uy*ic, uy*uz*ic - ux*s],
        [uz*ux*ic - uy*s, uz*uy*ic + ux*s, c + uz*uz*ic]
    ], dtype=float)
    pivot_np = np.array(pivot.get_array(), dtype=float)
    for res in structure:
        for atom in res.get_atoms():
            coord = np.array(atom.get_coord(), dtype=float)
            shifted = coord - pivot_np
            rotated = R @ shifted
            atom.set_coord(rotated + pivot_np)


def bend(dna, n_hinge, alpha, gamma, zeta, epsilon):
    if isinstance(dna, (list, tuple)):
        residues = dna
    else:
        model = dna[0] if hasattr(dna, '__getitem__') else dna
        chain = next(model.get_chains())
        residues = list(chain.get_residues())
    residues = rotate_alpha(residues, alpha, n_hinge)
    residues = rotate_gamma(residues, gamma, n_hinge)
    residues = rotate_zeta(residues, zeta, n_hinge)
    residues = rotate_epsilon(residues, epsilon, n_hinge)
    return residues


def bend_head(dna, n_hinge, alpha, gamma, zeta, epsilon):
    if isinstance(dna, (list, tuple)):
        residues = dna
    else:
        model = dna[0] if hasattr(dna, '__getitem__') else dna
        chain = next(model.get_chains())
        residues = list(chain.get_residues())
    residues = rotate_head_alpha(residues, alpha, n_hinge)
    residues = rotate_head_gamma(residues, gamma, n_hinge)
    residues = rotate_head_zeta(residues, zeta, n_hinge)
    residues = rotate_head_epsilon(residues, epsilon, n_hinge)
    return residues


def check_geometry(structure, res_per_seg, nedge):
    total_loss = 0.0
    for chain in structure[0]:
        residues = list(chain.get_residues())
        n = len(residues)
        if n == 0:
            continue
        for seg_end in range(res_per_seg - 1, n, res_per_seg):
            if seg_end + 1 >= n:
                continue
            total_loss += _single_joint_loss(residues[seg_end], residues[seg_end + 1])
        last_res = residues[-1]
        first_res = residues[0]
        total_loss += _single_joint_loss(last_res, first_res)
    return total_loss


def _single_joint_loss(res_up, res_down):
    try:
        o3_up  = res_up['O3\''].get_coord()
        c3_up  = res_up['C3\''].get_coord()
        p_down = res_down['P'].get_coord()
        o5_down= res_down['O5\''].get_coord()
    except KeyError:
        return 0.0
    distance = np.linalg.norm(p_down - o3_up)
    angle1   = angle_between_vectors(o3_up - p_down, o3_up - c3_up)
    angle2   = angle_between_vectors(p_down - o3_up, p_down - o5_down)
    loss = (10.0 * (distance - 1.6)**2 + 0.01 * (angle1 - 119.0)**2 + 0.01 * (angle2 - 101.4)**2)
    return loss


# ========== 4. 主构建函数 ==========

def build_tetrahedron(sequences: dict, L: float, edge_rots: list,
                      per_seg_bends: list, *, out_pdb: str = 'output.pdb',
                      len_hinge: int = 1):
    if len(edge_rots) != 6:
        raise ValueError('edge_rots 必须恰好 6 个浮点数')
    if len(per_seg_bends) != 12:
        raise ValueError('per_seg_bends 必须恰好 12 个元素')
    if any(len(row) != 4 for row in per_seg_bends):
        raise ValueError('每个 bend 元素必须是 4 个 float')

    tmp_pdb = '_tmp_12seg.pdb'
    seg_len = len(sequences['A']) // 3

    gen_segments_multi(list(sequences.items()), 3, tmp_pdb)

    parser = PDBParser(QUIET=True)
    st12 = parser.get_structure('z12', tmp_pdb)
    offset_z = L/2 - 3.38 * (seg_len - 1) / 2
    for atom in st12.get_atoms():
        x, y, z = atom.get_coord()
        atom.set_coord([x, y, z + offset_z])

    core = TetrahedronCore(sequences, ignore_last_bases=len_hinge)
    ass = core.validate_edges()
    
    
    ABSOLUTE_ANGLE_DEF = [
        ('F', 'EF', 'FG'), ('G', 'FG', 'GE'), ('E', 'GE', 'EF'),
        ('H', 'EH', 'FH'), ('H', 'FH', 'GH'), ('H', 'GH', 'EH'),
        ('E', 'EF', 'EH'), ('E', 'EH', 'GE'), ('F', 'FG', 'FH'),
        ('F', 'FH', 'EF'), ('G', 'GE', 'GH'), ('G', 'GH', 'FG'),
    ]
    
    seg_to_edge = {}
    for edge_name, info in ass.items():
        ch1, ch2 = info['chains']
        seg1, seg2 = info['segments']
        seg_to_edge[(ch1, seg1)] = edge_name
        seg_to_edge[(ch2, seg2)] = edge_name
    
    angle_to_dna = {}
    for angle_idx, (vertex, edge_in, edge_out) in enumerate(ABSOLUTE_ANGLE_DEF):
        found = False
        for ch in 'ABCD':
            for seg_idx in range(3):
                edge = seg_to_edge.get((ch, seg_idx))
                if edge is None:
                    continue
                if edge == edge_in:
                    next_seg = (seg_idx + 1) % 3
                    next_edge = seg_to_edge.get((ch, next_seg))
                    if next_edge == edge_out:
                        angle_to_dna[angle_idx] = [
                            (ch, seg_idx, 'tail'),
                            (ch, next_seg, 'head')
                        ]
                        found = True
                        break
            if found:
                break
    
    per_seg_bends_map = {}
    for angle_idx, positions in angle_to_dna.items():
        for ch, seg_idx, position in positions:
            key = (ch, seg_idx)
            if key not in per_seg_bends_map:
                per_seg_bends_map[key] = [None, None]
            if position == 'head':
                per_seg_bends_map[key][0] = per_seg_bends[angle_idx]
            else:
                per_seg_bends_map[key][1] = per_seg_bends[angle_idx]

    def _dock_and_twist(structure, ass, seg_len, len_hinge, edge_rots):
        FIXED_EDGE_ORDER = ['EF', 'FG', 'GE', 'EH', 'FH', 'GH']
        model = structure[0]
        chain_res = {ch.id: list(ch.get_residues()) for ch in model}
        nres_full = seg_len
        z_axis = Vector(0, 0, 1)
        if len_hinge == 1 or len_hinge == 3:
            dz = 3.38 * 1###这里记得修改
        if len_hinge == 2 or len_hinge == 4:
            dz = 0
        dphi = -85/2 + 36*(seg_len-3)###这里记得修改
        #print(dphi-360*(dphi//360))
        #print(seg_len)
        comp_map = {}
        for edge, info in ass.items():
            c1, c2 = info['chains']
            i1, i2 = info['segments']
            if c1 < c2:
                comp_map[(c2, i2)] = (c1, i1)
            else:
                comp_map[(c1, i1)] = (c2, i2)
        
        for (t_ch, t_idx), (s_ch, s_idx) in comp_map.items():
            mobile = chain_res[s_ch if t_ch < s_ch else t_ch][
                     (s_idx if t_ch < s_ch else t_idx)*nres_full :
                     (s_idx if t_ch < s_ch else t_idx+1)*nres_full]
            for res in mobile:
                for atom in res:
                    atom.transform(np.eye(3), np.array([0, 0, dz]))
            rotate_around_axis(mobile, -dphi, z_axis, Vector(0,0,0))

        for idx, edge in enumerate(FIXED_EDGE_ORDER):
            info = ass[edge]
            c1, c2 = info['chains']
            i1, i2 = info['segments']
            seg1 = chain_res[c1][i1*nres_full:(i1+1)*nres_full]
            seg2 = chain_res[c2][i2*nres_full:(i2+1)*nres_full]
            rotate_around_axis(seg1, edge_rots[idx], z_axis, Vector(0,0,0))
            rotate_around_axis(seg2, -edge_rots[idx], z_axis, Vector(0,0,0))
    
    _dock_and_twist(st12, ass, seg_len, len_hinge, edge_rots)


    chain_order = [ch for ch in 'ABCD' for _ in range(3)]
    
    # 根据 len_hinge 确定弯折参数
    if len_hinge == 1 or len_hinge == 3:
        n_hinge_head = 1
        n_hinge_tail = 2
    elif len_hinge == 2 or len_hinge == 4:
        n_hinge_head = 2
        n_hinge_tail = 2
    else:
        raise ValueError(f"不支持的 len_hinge={len_hinge}，仅支持 1-4")

    # 尾部弯折（3'端）
    for seg_idx, ch_id in enumerate(chain_order):
        chain = next(c for c in st12[0] if c.id == ch_id)
        residues = list(chain.get_residues())
        start = seg_idx % 3 * seg_len
        end = start + seg_len
        sub_res = residues[start:end]
        angles = per_seg_bends_map[(ch_id, seg_idx % 3)][1]
        if n_hinge_tail > 0:  # len_hinge=3时执行，len_hinge=1时也执行
            sub_res = bend(sub_res, n_hinge_tail, alpha=angles[0], gamma=angles[1],
                           zeta=angles[2], epsilon=angles[3])
    
    # 头部弯折（5'端）
    for seg_idx, ch_id in enumerate(chain_order):
        chain = next(c for c in st12[0] if c.id == ch_id)
        residues = list(chain.get_residues())
        start = seg_idx % 3 * seg_len
        end = start + seg_len
        sub_res = residues[start:end]
        angles = per_seg_bends_map[(ch_id, seg_idx % 3)][0]
        if n_hinge_head > 0:  # len_hinge=1时执行，len_hinge=3时跳过
            sub_res = bend_head(sub_res, n_hinge_head, alpha=angles[0], gamma=angles[1],
                                zeta=angles[2], epsilon=angles[3])
    
    center = Vector(np.sqrt(3)/6 * L, 0., L/2.)
    y_axis = Vector(0, 1, 0)

    chain_match_idx = {}
    A_base_idx = {}
    for edge, info in ass.items():
        a_ch, a_idx = info['chains'][0], info['segments'][0]
        o_ch, o_idx = info['chains'][1], info['segments'][1]
        if a_ch == 'A':
            chain_match_idx[o_ch] = o_idx
            A_base_idx[o_ch] = a_idx

    for ch_id in 'ABCD':
        chain = next(c for c in st12[0] if c.id == ch_id)
        Cn_rotate(chain, 3, y_axis, center)
        if ch_id != 'A' and ch_id in chain_match_idx:
            offset = chain_match_idx[ch_id] - A_base_idx[ch_id]
            if offset != 0:
                angle = offset * (-120.0)
                rotate_around_axis(chain, angle, y_axis, center)

    a_axis_formula = [
        (Vector(0, 0, 0), Vector(0, 0, L)),
        (Vector(0, 0, L), Vector(np.sqrt(3)*L/2, 0, L/2)),
        (Vector(np.sqrt(3)*L/2, 0, L/2), Vector(0, 0, 0))
    ]
    non_a_axis = {}
    for edge, info in ass.items():
        a_ch, a_idx = info['chains'][0], info['segments'][0]
        o_ch, o_idx = info['chains'][1], info['segments'][1]
        if o_ch == 'A': continue
        non_a_axis.setdefault(o_ch, a_idx)

    for ch_id in 'BCD':
        if ch_id not in non_a_axis: continue
        A = Vector(np.sqrt(3)*L/6, np.sqrt(6)*L/3, L/2)
        a_idx = non_a_axis[ch_id]
        start, end = a_axis_formula[a_idx]
        axis_vec = (end - start).normalized()
        axis_vec2 = (A - (start+end)/2).normalized()
        chain = next(c for c in st12[0] if c.id == ch_id)
        sele = list(chain.get_residues())
        rotate_around_axis(sele, 70.53, axis_vec, start)
        rotate_around_axis(sele, 180., axis_vec2, A)

    io = PDBIO()
    io.set_structure(st12)
    total_loss = check_geometry(st12, res_per_seg=seg_len, nedge=3)
    
    if out_pdb is not None:
        io.save(out_pdb)
    
    if os.path.exists(tmp_pdb):
        os.remove(tmp_pdb)
    #print("total_loss="+str(total_loss))
    return total_loss

if __name__ == "__main__":
    import sys
    import time
    import os

    # 解析命令行参数
    args = sys.argv[1:]
    use_autofill = '-autofill' in args

    # 去掉 -autofill 参数，获取文件路径
    file_args = [a for a in args if a != '-autofill']

    if not file_args:
        print("用法:")
        print("  标准模式:  python TDN_builder.py <序列文件.txt>")
        print("  自动补全:  python TDN_builder.py -autofill <序列文件.txt>")
        print("")
        print("-autofill 模式输入格式:")
        print("  第1行: len_hinge = N")
        print("  第2行: 完整A链（3段）")
        print("  第3行: B[0]+B[1]（2段）")
        print("  第4行: C[0]（1段）")
        print("  程序自动补全: B[2]=revc(A[0]), C[1]=revc(A[1]), C[2]=revc(B[1])")
        print("                D[0]=revc(A[2]), D[1]=revc(C[0]), D[2]=revc(B[0])")
        sys.exit(1)

    txt_file = file_args[0]

    print(f"Reading from: {txt_file}")

    # 如果启用 -autofill 模式
    if use_autofill:
        print("\n[Autofill Mode] 自动补全序列...")

        # 检测是否误用了标准格式
        is_std, msg = check_input_mismatch(txt_file)
        if is_std:
            print(f"\n警告: {msg}")
            print("\n如果您确认要使用 -autofill 模式，请按 Enter 继续...")
            print("或按 Ctrl+C 取消并修改输入文件。")
            try:
                input()
            except KeyboardInterrupt:
                print("\n已取消。")
                sys.exit(1)

        # 处理自动补全
        try:
            standard_file = process_autofill(txt_file, out_path=None, verbose=True)
            txt_file = standard_file  # 使用补全后的文件继续处理
            print(f"\n使用补全后的文件继续: {txt_file}")
        except ValueError as e:
            print(f"\n[Autofill Error] {e}")
            sys.exit(1)

    # 1. 读取配置和序列（原始链或补全后的链）
    len_hinge, seq_a, seq_2, seq_3, seq_4 = read_config_and_sequences(txt_file)

    raw_seqs = {
        'A': seq_a,
        'B': seq_2,
        'C': seq_3,
        'D': seq_4
    }
    result = find_valid_configuration(raw_seqs, len_hinge)
    
    if result is None:
        seqs = auto_assign_sequences(seq_a, [seq_2, seq_3, seq_4], 
                                     ignore_last_bases=len_hinge)
        cut_positions = None
    else:
        partition, corrected_seqs, cut_positions = result
        #print(f"[V2] 段分配: {partition['assignment']}")
        seqs = corrected_seqs  # 使用修正链构建DNA
    
    # 3. 验证拓扑
    try:
        from pair_DNA_tetrahedron_v2 import TetrahedronCore
        core = TetrahedronCore(seqs, ignore_last_bases=len_hinge, verbose=False)
        ass = core.validate_edges()
        print(f"\nsequences paired! (hinge={len_hinge})")
    except Exception as e:
        print(f"\n拓扑验证失败: {e}")
        sys.exit(1)
    
    # 4. 运行构建
    a = time.time()
    
    L, edge_rots, per_seg_bends = dihedrals_and_axis(seqs['A'], len_hinge)
    
    tmp_pdb = 'tetrahedron_tmp.pdb'
    final_pdb = 'tetrahedron_result.pdb'
    
    total_loss = build_tetrahedron(
        seqs, L, edge_rots, per_seg_bends,
        out_pdb=tmp_pdb,
        len_hinge=len_hinge
    )
    
    # 5. 调用rearrangeDNA调回原始顺序（如果V2优化成功）
    if cut_positions:
        for ch in ['A', 'B', 'C', 'D']:
            seq = seqs[ch]
            cut_pos = cut_positions[ch]
        try:
            # 导入 rearrangeDNA 中的函数
            import importlib.util
    
            # 加载 rearrangeDNA 模块
            spec = importlib.util.spec_from_file_location("rearrangeDNA_v2", "rearrangeDNA_v2.py")
            rearrange_module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(rearrange_module)
    
            # 使用 rearrange_all_chains 对所有链进行重排
            success = rearrange_module.rearrange_all_chains(tmp_pdb, final_pdb, cut_positions)
            
            if success:
                #print(f"\n[V2] 调整完成，已调回原始顺序: {final_pdb}")
                if os.path.exists(tmp_pdb):
                    os.remove(tmp_pdb)
            else:
                print(f"[V2] 调整失败，使用修正链结构")
                os.rename(tmp_pdb, final_pdb)
        except Exception as e:
            print(f"[V2] rearrangeDNA调用失败: {e}")
            os.rename(tmp_pdb, final_pdb)
    else:
        os.rename(tmp_pdb, final_pdb)
    strip_terminal_phosphate(final_pdb)
    fix_ter_placement(final_pdb)
    b = time.time()
    
    # 6. 输出最终信息
    #print(f"\n{'='*60}")
    #print(f"[V2] 最终输出:")
    #if cut_positions:
        #print(f"  结构已调回原始链顺序")
        #print(f"  剪切位置: {cut_positions}")
    #else:
        #print(f"  使用原始配置（未优化）")
    #print(f"{'='*60}")
    
    print(f"\ntotal_loss={total_loss}")
    print(f"total time: {round(b-a, 3)}s")
    print(f"output file: {final_pdb}")
