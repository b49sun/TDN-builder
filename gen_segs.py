# -*- coding: utf-8 -*-
"""
build_63mer_from_3x21.py
生成 3 条 21-mer 直链 → 残基号 1–63 → 同一条 Chain('A')
不旋转、不平移，坐标与单条 21-mer 完全一致
"""
import tempfile, os
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain
from generate_dna import gen_DNA   # 你的接口


def gen_segments_multi(sequences_items, n_edge, out_pdb):
    """
    sequences_items : list[tuple[str, str]]  [(链名, 序列), ...] 如 [('A', 'ATCG...'), ('B', '...'), ...]
    n_edge          : int        每条链要切成的段数
    out_pdb         : str        输出 PDB 文件名
    """
    if n_edge <= 0:
        raise ValueError("n_edge 必须是正整数")

    if not sequences_items:
        raise ValueError("请至少提供一条序列")

    # 第一条序列的校验
    first_name, first_seq = sequences_items[0]
    first_len = len(first_seq)
    if first_len % n_edge != 0:
        raise ValueError(f"第一条序列长度 {first_len} 不能被 n_edge={n_edge} 整除")
    len_per_edge = first_len // n_edge

    # 其余序列长度校验
    for idx, (name, seq) in enumerate(sequences_items[1:], start=2):
        if len(seq) != first_len:
            raise ValueError(f"第 {idx} 条序列（{name}）长度 {len(seq)} 与第一条 {first_len} 不等")

    # 构建 PDB 框架
    structure = Structure.Structure("ssDNA")
    model = Model.Model(0)
    structure.add(model)

    # 为每条序列生成独立 Chain
    for chain_id, seq in sequences_items:  # 直接使用传入的链名
        chain = Chain.Chain(chain_id)       # 不再是 chr(ord('A')+idx)
        model.add(chain)

        # 生成 n_edge 段完全重合的片段
        for seg_idx in range(n_edge):
            tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
            tmp_name = tmp.name
            tmp.close()
            start   = seg_idx * len_per_edge
            end     = (seg_idx + 1) * len_per_edge
            seg_seq = seq[start:end]
            # 生成单段直线（坐标永远相同）
            gen_DNA(tmp_name, seg_seq, 36.0, 3.38, single=True)

            # 读入
            parser = PDBParser(QUIET=True)
            st = parser.get_structure('tmp', tmp_name)
            residues = list(next(st[0].get_chains()).get_residues())

            # 残基编号顺次接龙
            base_num = seg_idx * len_per_edge + 1
            for res in residues:
                res.id = (res.id[0], base_num, res.id[2])
                chain.add(res)
                base_num += 1

            os.remove(tmp_name)

    # 写出 PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_pdb)


# ------------------ CLI 示例 ------------------
if __name__ == '__main__':
    # 示例输入，现在需要显式指定链名
    sequences_items = [
        ('A', 63 * 'A'),
        ('B', 63 * 'T'),
    ]
    n_edge = 3
    out_file = 'multi_chain_DNA.pdb'
    gen_segments_multi(sequences_items, n_edge, out_file)