"""
DNA序列随机生成器
生成3条ATCG序列，长度比为3:2:1，长度和AT比例可自选

使用方法:
    from dna_generator import generate_dna_sequences

    # 示例：基准长度200，AT比例60%
    seqs, stats = generate_dna_sequences(base_length=200, at_ratio=0.6)

    for i, seq in enumerate(seqs):
        print(f"序列{i+1} (长度 {len(seq)}): {seq}")
"""

import random


def generate_dna_sequences(length_ratio_3=3, length_ratio_2=2, length_ratio_1=1, 
                           base_length=100, at_ratio=0.5, seed=None):
    """
    随机生成3条DNA序列，长度比为3:2:1

    参数:
        length_ratio_3: 第一条序列的长度比例（默认3）
        length_ratio_2: 第二条序列的长度比例（默认2）
        length_ratio_1: 第三条序列的长度比例（默认1）
        base_length: 基准长度（比例1对应的实际长度，默认100）
        at_ratio: AT比例（0-1之间，默认0.5表示AT和CG各占50%）
        seed: 随机种子（可选，用于复现结果）

    返回:
        sequences: 包含3条序列的列表 [seq1, seq2, seq3]
        stats: 各序列的详细统计信息列表
    """
    if seed is not None:
        random.seed(seed)

    # 计算各序列长度
    len1 = base_length * length_ratio_3
    len2 = base_length * length_ratio_2
    len3 = base_length * length_ratio_1

    # 根据AT比例确定各碱基的概率
    # A和T各占 at_ratio/2，C和G各占 (1-at_ratio)/2
    p_a = at_ratio / 2
    p_t = at_ratio / 2
    p_c = (1 - at_ratio) / 2
    p_g = (1 - at_ratio) / 2

    bases = ['A', 'T', 'C', 'G']
    weights = [p_a, p_t, p_c, p_g]

    # 生成序列
    seq1 = ''.join(random.choices(bases, weights=weights, k=len1))
    seq2 = ''.join(random.choices(bases, weights=weights, k=len2))
    seq3 = ''.join(random.choices(bases, weights=weights, k=len3))

    sequences = [seq1, seq2, seq3]

    # 统计各序列的碱基组成
    stats = []
    for i, seq in enumerate(sequences):
        a_count = seq.count('A')
        t_count = seq.count('T')
        c_count = seq.count('C')
        g_count = seq.count('G')
        at_count = a_count + t_count
        cg_count = c_count + g_count
        at_pct = at_count / len(seq) * 100 if len(seq) > 0 else 0

        stats.append({
            'seq_id': i + 1,
            'length': len(seq),
            'A': a_count,
            'T': t_count,
            'C': c_count,
            'G': g_count,
            'AT_count': at_count,
            'CG_count': cg_count,
            'AT_ratio': f"{at_pct:.1f}%"
        })

    return sequences, stats


def print_sequences(sequences, stats):
    """美观地打印序列和统计信息"""
    print("=" * 60)
    print("DNA序列生成结果")
    print("=" * 60)
    for s, st in zip(sequences, stats):
        print(f"\n序列{st['seq_id']} (长度: {st['length']}):")
        # 每行显示60个碱基，方便阅读
        for j in range(0, len(s), 60):
            print(f"  {s[j:j+60]}")
        print(f"  A:{st['A']} T:{st['T']} C:{st['C']} G:{st['G']} "
              f"| AT:{st['AT_count']} CG:{st['CG_count']} | AT占比: {st['AT_ratio']}")
    print("=" * 60)


# ============ 主程序入口 ============
if __name__ == "__main__":
    seqs, stats = generate_dna_sequences(base_length=19, at_ratio=0.5, seed=42)
    print_sequences(seqs, stats)
