# -*- coding: utf-8 -*-
from Bio.PDB.vectors import rotaxis, Vector
from Bio.PDB import PDBParser, PDBIO
from generate_dna import gen_DNA
from gen_segs import gen_segments
from local_math import angle_between_vectors, dihedral_among_vectors, get_normal_vector
import numpy as np
import math
import time

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# triangle_DNA.py

# ---------- 参数 ----------

# 4. 分段旋转

def Cn_rotate(structure,
                              n: int,
                              axis: Vector,
                              p0: Vector = Vector(0, 0, 0)) -> None:
    if len(structure) % n:
        raise ValueError('残基数必须能被 n 整除')
    total_res = len(structure)
    res_per_seg = total_res // n

    # ---- 工具：Rodrigues 矩阵 ----
    def rot_mat(unit_ax: np.ndarray, theta: float) -> np.ndarray:
        c, s = np.cos(theta), np.sin(theta)
        ic = 1.0 - c
        ux, uy, uz = unit_ax
        return np.array([
            [c + ux*ux*ic,   ux*uy*ic - uz*s, ux*uz*ic + uy*s],
            [uy*ux*ic + uz*s, c + uy*uy*ic,   uy*uz*ic - ux*s],
            [uz*ux*ic - uy*s, uz*uy*ic + ux*s, c + uz*uz*ic]
        ], dtype=float)

    z_axis = np.array([0., 0., 1.])
    axis_np = np.array(axis.get_array(), dtype=float)
    axis_unit = axis_np / np.linalg.norm(axis_np)
    p0_np = np.array(p0.get_array(), dtype=float)

    # ========== 3. 绕指定轴 Cn 旋转 ==========
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

def rotate_alpha(structure, nedge, theta_deg):
    theta = -math.radians(theta_deg)
    model = structure[0]
    chain = next(model.get_chains())
    residues = list(chain.get_residues())
    nres_in_edge = len(residues) // nedge

    for i in range(nedge):
        bend_at = nres_in_edge * (i + 1)-2
        # 旋转轴：C5' -> C4'
        atom_P = residues[bend_at]['P']
        atom_O5 = residues[bend_at]['O5\'']
        axis_vec = Vector(atom_O5.get_coord()) - Vector(atom_P.get_coord())
        pivot = Vector(atom_P.get_coord())

        # 旋转矩阵
        rot_mat = rotaxis(theta, axis_vec)

        # 只转当前残基的 C4' 及之后原子
        for atom in residues[bend_at+1]:
            #if atom.get_name() in ['C5\'', 'O5\'', 'P', 'OP1', 'OP2']:
            coord = np.array(atom.get_coord())
            translated = coord - pivot.get_array()
            rotated = rot_mat @ translated
            new_coord = rotated + pivot.get_array()
            atom.set_coord(new_coord)
        for atom in residues[bend_at]:
            if atom.get_name() in ['OP1', 'OP2']:
                coord = np.array(atom.get_coord())
                translated = coord - pivot.get_array()
                rotated = rot_mat @ translated
                new_coord = rotated + pivot.get_array()
                atom.set_coord(new_coord)
    return structure

def rotate_gamma(structure, nedge, theta_deg):#epsilon
    theta = -math.radians(theta_deg)
    model = structure[0]
    chain = next(model.get_chains())
    residues = list(chain.get_residues())
    nres_in_edge = len(residues) // nedge

    for i in range(nedge):
        bend_at = nres_in_edge * (i + 1)-1
        # 旋转轴：C5' -> C4'
        atom_C5 = residues[bend_at]['C5\'']
        atom_C4 = residues[bend_at]['C4\'']
        axis_vec = Vector(atom_C4.get_coord()) - Vector(atom_C5.get_coord())
        pivot = Vector(atom_C5.get_coord())

        # 旋转矩阵
        rot_mat = rotaxis(theta, axis_vec)

        # 只转当前残基的 C4' 及之后原子
        for atom in residues[bend_at]:
            if atom.get_name() in ['C5\'', 'O5\'', 'P', 'OP1', 'OP2']:
                coord = np.array(atom.get_coord())
                translated = coord - pivot.get_array()
                rotated = rot_mat @ translated
                new_coord = rotated + pivot.get_array()
                atom.set_coord(new_coord)
    return structure

def rotate_zeta(structure, nedge, theta_deg):
    theta = -math.radians(theta_deg)
    model = structure[0]
    chain = next(model.get_chains())
    residues = list(chain.get_residues())
    nres_in_edge = len(residues) // nedge

    for i in range(nedge):
        bend_at = nres_in_edge * (i + 1)-1
        # 旋转轴：C5' -> C4'
        atom_O3 = residues[bend_at]['O3\'']
        atom_P = residues[bend_at-1]['P']
        axis_vec = Vector(atom_P.get_coord())-Vector(atom_O3.get_coord())
        pivot = Vector(atom_O3.get_coord())

        # 旋转矩阵
        rot_mat = rotaxis(theta, axis_vec)

        # 只转当前残基的 C4' 及之后原子
        for atom in residues[bend_at]:
            if atom.get_name() not in ['O3\'']:
                coord = np.array(atom.get_coord())
                translated = coord - pivot.get_array()
                rotated = rot_mat @ translated
                new_coord = rotated + pivot.get_array()
                atom.set_coord(new_coord)
    
    return structure

def rotate_epsilon(structure, nedge, theta_deg):
    theta = -math.radians(theta_deg)
    model = structure[0]
    chain = next(model.get_chains())
    residues = list(chain.get_residues())
    nres_in_edge = len(residues) // nedge

    for i in range(nedge):
        bend_at = nres_in_edge * (i + 1)-1
        # 旋转轴：C5' -> C4'
        atom_C3 = residues[bend_at]['C3\'']
        atom_O3 = residues[bend_at]['O3\'']
        axis_vec = Vector(atom_O3.get_coord()) - Vector(atom_C3.get_coord())
        pivot = Vector(atom_C3.get_coord())

        # 旋转矩阵
        rot_mat = rotaxis(theta, axis_vec)

        # 只转当前残基的 C4' 及之后原子
        for atom in residues[bend_at]:
            if atom.get_name() not in ['O3\'']:
                coord = np.array(atom.get_coord())
                translated = coord - pivot.get_array()
                rotated = rot_mat @ translated
                new_coord = rotated + pivot.get_array()
                atom.set_coord(new_coord)
    return structure

def rotate_around_axis(structure, angle_deg: float,
                       axis_vec, pivot) -> None:
    """
    将 structure 中所有原子绕 pivot 点、沿 axis_vec 方向旋转 angle_deg 度。
    右手定则：拇指指向 axis_vec 正向，四指弯曲方向为正角度方向。
    """
    # 自动转成 Vector
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

    return structure

def bend(dna, nedge, alpha, gamma, zeta, epsilon):
    dna = rotate_alpha(dna, nedge, alpha)
    dna = rotate_gamma(dna, nedge, gamma)
    dna = rotate_zeta(dna, nedge, zeta)
    dna = rotate_epsilon(dna, nedge, epsilon)
    return dna

'''
def check_geometry(structure,res_per_seg):
    # 返回折叠后的首位原子距离

    model = structure[0]
    chain = next(model.get_chains())
    chain.id = 'A'
    residues = list(chain.get_residues())
    ref_O3 = residues[res_per_seg]['O3\''].get_coord()
    ref_C3 = residues[res_per_seg]['C3\''].get_coord()
    ref_C4 = residues[res_per_seg]['C4\''].get_coord()
    atom_P = residues[res_per_seg-1]['P'].get_coord()
    atom_O5 = residues[res_per_seg-1]['O5\''].get_coord()
    atom_C5 = residues[res_per_seg-1]['C5\''].get_coord()
    distance_P_O3 = np.linalg.norm(atom_P - ref_O3)
    angle_C3_O3_P = angle_between_vectors(ref_O3-atom_P,ref_O3-ref_C3)
    dihedral_C4_C3_O3_P = dihedral_among_vectors(ref_O3-atom_P,ref_C3-ref_O3,ref_C4-ref_C3)
    angle_O3_P_O5 = angle_between_vectors(atom_P-ref_O3,atom_P-atom_O5)
    dihedral_C3_O3_P_O5 = dihedral_among_vectors(atom_P-atom_O5,ref_O3-atom_P,ref_C3-ref_O3)
    dihedral_O3_P_O5_C5 = dihedral_among_vectors(atom_O5-atom_C5,atom_P-atom_O5,ref_O3-atom_P)
    loss = 10*(distance_P_O3-1.6)**2 + 0.01*(angle_C3_O3_P-119)**2+0.01*(angle_O3_P_O5-101.4)**2
    #+0.01*(dihedral_C4_C3_O3_P-159.1)**2+0.01*(dihedral_C3_O3_P_O5+98.9)**2+0.01*(dihedral_O3_P_O5_C5+39.2)**2
    return loss
'''
def check_geometry(structure, res_per_seg, nedge):
    """
    返回所有转角损失之和
    """
    from local_math import angle_between_vectors   # 确保导入
    model = structure[0]
    chain = next(model.get_chains())
    residues = list(chain.get_residues())

    total_loss = 0.0
    for idx in range(1, nedge + 1):
        center_idx = idx * res_per_seg
        if center_idx >= len(residues):
            break

        ref_O3  = residues[center_idx]['O3\''].get_coord()
        ref_C3  = residues[center_idx]['C3\''].get_coord()
        atom_P  = residues[center_idx - 1]['P'].get_coord()
        atom_O5 = residues[center_idx - 1]['O5\''].get_coord()

        distance_P_O3 = np.linalg.norm(atom_P - ref_O3)
        angle_C3_O3_P = angle_between_vectors(ref_O3 - atom_P, ref_O3 - ref_C3)
        angle_O3_P_O5 = angle_between_vectors(atom_P - ref_O3, atom_P - atom_O5)

        loss = (10.0 * (distance_P_O3 - 1.6) ** 2 +
                0.01 * (angle_C3_O3_P - 119.0) ** 2 +
                0.01 * (angle_O3_P_O5 - 101.4) ** 2)
        total_loss += loss

    return total_loss

if __name__ == "__main__":
    a = time.time()
    nedge, n_in_edge = 3,3
    seq_len = 10 * n_in_edge +1         # 总碱基数
    input_pdb  = 'single_DNA.pdb'
    output_pdb = 'triangle_DNA.pdb'

    gen_segments(seq_len, nedge, input_pdb)

# 2. 读取结构
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('ssDNA', input_pdb)
    print(structure[0][0])
    model = structure[0]
    chain = list(model.get_chains())[0]   # 单链
    for atom in structure.get_atoms():
        atom.occupancy = 1.0
# 3. 把原子按残基序号排序，确保顺序正确
    all_residues = sorted(chain.get_residues(), key=lambda r: r.id[1])
    total_res = len(all_residues)
    res_per_seg = total_res // nedge        # 每段残基数
    if seq_len == 10 * n_in_edge+1 and nedge ==3:
        bend(structure,nedge,-1.2458673512666154,-26.463643309290223,-24.47749822011577,179.99993143030852)
        Cn_rotate(all_residues, nedge, Vector(0, 1, 0),Vector(21.47+(n_in_edge-2)*33.8/(2*math.sqrt(3)), 0.00, 36.60+(n_in_edge-2)*33.8/2))
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        print(f"Triangle DNA written to {output_pdb}")
        print("loss = "+str(check_geometry(structure,seq_len,nedge)))
    if seq_len == 10 * n_in_edge and nedge ==3:
        bend(structure,nedge,-76.24092011815492, -23.04459312769806, -129.99394894329447, 33.44059512498262)
        Cn_rotate(all_residues, nedge, Vector(0, 1, 0),Vector(24.97+(n_in_edge-2)*33.8/(2*math.sqrt(3)), 0.00, 33.93+(n_in_edge-2)*33.8/2))
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        print(f"Triangle DNA written to {output_pdb}")
        print("loss = "+str(check_geometry(structure,seq_len,nedge)))
    if seq_len == 10 * n_in_edge+1 and nedge ==4:
        bend(structure,nedge,-79.99277107791204, -22.10650792765256, -125.33661887849328, 0.45587744954605985)
        Cn_rotate(all_residues, nedge, Vector(0, 1, 0),Vector(38+(n_in_edge-2)*33.8/2, 0.00, 37.48+(n_in_edge-2)*33.8/2))
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        print(f"Triangle DNA written to {output_pdb}")
        print("loss = "+str(check_geometry(structure,seq_len,nedge)))
    if seq_len == 10 * n_in_edge and nedge ==4:
        bend(structure,nedge,-0.7082904508924249, -21.362550689180946, -0.4209239908627711, -5.0)
        Cn_rotate(all_residues, nedge, Vector(0, 1, 0),Vector(38+(n_in_edge-2)*33.8/2, 0.00, 34.17+(n_in_edge-2)*33.8/2))
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        print(f"Triangle DNA written to {output_pdb}")
        print("loss = "+str(check_geometry(structure,seq_len,nedge)))
    if seq_len == 10 * n_in_edge+1 and nedge ==5:
        bend(structure,nedge,-5.045250974712101, -15.032755225433865, 13.258953672229154, 39.13752458909933)
        Cn_rotate(all_residues, nedge, Vector(0, 1, 0),Vector(48.26+(n_in_edge-2)*33.8*1.3764/2, 0.00, 39.17+(n_in_edge-2)*33.8/2))
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        print(f"Triangle DNA written to {output_pdb}")
        print("loss = "+str(check_geometry(structure,seq_len,nedge)))
    if seq_len == 10 * n_in_edge and nedge ==5:
        bend(structure,nedge,-26.39215412321299, 5.971944357361006, 17.049770870774832, 6.398366655573603)
        Cn_rotate(all_residues, nedge, Vector(0, 1, 0),Vector(46.92+(n_in_edge-2)*33.8*1.3764/2, 0.00, 33.88+(n_in_edge-2)*33.8/2))
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        print(f"Triangle DNA written to {output_pdb}")
        print("loss = "+str(check_geometry(structure,seq_len,nedge)))
    if seq_len == 10 * n_in_edge and nedge ==6:
        bend(structure,nedge,14.506596795505912, 26.694220994952317, -16.802790680041937, 27.196029202212184)
        Cn_rotate(all_residues, nedge, Vector(0, 1, 0),Vector(62.71+(n_in_edge-2)*33.8*np.sqrt(3)/2, 0.00, 34.12+(n_in_edge-2)*33.8/2))
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        print(f"Triangle DNA written to {output_pdb}")
        print("loss = "+str(check_geometry(structure,seq_len,nedge)))
    if seq_len == 10 * n_in_edge+1 and nedge ==6:
        bend(structure,nedge,-53.56404654583809, -29.086221426542927, -127.05802482114856, 49.20422551111347)
        Cn_rotate(all_residues, nedge, Vector(0, 1, 0),Vector(66.45+(n_in_edge-2)*33.8*np.sqrt(3)/2, 0.00, 37.91+(n_in_edge-2)*33.8/2))
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        print(f"Triangle DNA written to {output_pdb}")
        print("loss = "+str(check_geometry(structure,seq_len,nedge)))
    b = time.time()
    print("生成结构耗时"+f'{(b-a):.2f}'+"s")