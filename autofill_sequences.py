# -*- coding: utf-8 -*-
"""
DNA 四面体序列自动补全模块 (Autofill Mode) - 纯文本处理，带 hinge 处理

输入格式（-autofill 模式）：
    第1行: len_hinge = N  (或单独一个数字N)
    第2行: 完整A链（3段）
    第3行: B[0] + B[1]（2段，B[2]自动生成）
    第4行: C[0]（1段，C[1]和C[2]自动生成）

hinge 处理规则：
    - 互补匹配时：忽略每段末尾 len_hinge 个碱基
    - 生成互补链：在末尾补 len_hinge 个 'A'

补全规则（6条边互补，忽略 hinge 区）：
    B[2] = revc(A[0][:-hinge]) + 'A'*hinge        # B[2]与A[0]互补（忽略hinge）
    C[1] = revc(A[1][:-hinge]) + 'A'*hinge        # C[1]与A[1]互补（忽略hinge）
    C[2] = revc(B[1][:-hinge]) + 'A'*hinge        # C[2]与B[1]互补（忽略hinge）
    D[0] = revc(A[2][:-hinge]) + 'A'*hinge        # D[0]与A[2]互补（忽略hinge）
    D[1] = revc(C[0][:-hinge]) + 'A'*hinge        # D[1]与C[0]互补（忽略hinge）
    D[2] = revc(B[0][:-hinge]) + 'A'*hinge        # D[2]与B[0]互补（忽略hinge）

输出: 补全为完整4条链的标准格式文件，之后与不开autofill一样处理
"""

import os
import sys
import re


def revc(seq: str) -> str:
    """反向互补: A<->T, C<->G"""
    comp = str.maketrans('ATCGatcg', 'TAGCtagc')
    return seq.translate(comp)[::-1]


def parse_autofill_input(txt_path: str):
    """
    读取-autofill模式的输入文件
    返回: (len_hinge, seq_a_full, seq_b_partial, seq_c_partial, seg_len)
    """
    with open(txt_path, 'r') as f:
        raw_lines = [line.strip() for line in f]

    if not raw_lines:
        raise ValueError("空文件")

    # 解析第一行：len_hinge
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

    # 收集序列（支持多行）
    sequence_data = raw_lines[1:]
    sequences_blocks = []
    current_block = []

    for line in sequence_data:
        if line == '':
            if current_block:
                sequences_blocks.append(''.join(current_block).upper())
                current_block = []
                if len(sequences_blocks) >= 3:
                    break
        else:
            current_block.append(line)

    if current_block and len(sequences_blocks) < 3:
        sequences_blocks.append(''.join(current_block).upper())

    # 如果只有3个非空行，直接作为3条链
    if len(sequences_blocks) != 3:
        non_empty = [l.upper() for l in sequence_data if l]
        if len(non_empty) == 3:
            sequences_blocks = non_empty

    if len(sequences_blocks) != 3:
        raise ValueError(
            f"-autofill模式需要3条序列（A完整+B2段+C1段），"
            f"解析到{len(sequences_blocks)}条"
        )

    seq_a = sequences_blocks[0]
    seq_b_input = sequences_blocks[1]
    seq_c_input = sequences_blocks[2]

    L = len(seq_a)
    if L % 3 != 0:
        raise ValueError(f"A链长度{L}不是3的倍数")

    seg_len = L // 3

    # 验证B输入长度 = 2段
    if len(seq_b_input) != 2 * seg_len:
        raise ValueError(
            f"B链输入长度应为{2*seg_len}（2段），实际为{len(seq_b_input)}"
        )

    # 验证C输入长度 = 1段
    if len(seq_c_input) != seg_len:
        raise ValueError(
            f"C链输入长度应为{seg_len}（1段），实际为{len(seq_c_input)}"
        )

    return len_hinge, seq_a, seq_b_input, seq_c_input, seg_len


def autofill_sequences(len_hinge: int, seq_a: str, seq_b_input: str, seq_c_input: str, seg_len: int):
    """
    根据输入序列自动补全所有4条链（带 hinge 处理）

    hinge 处理：
        - 互补匹配时忽略每段末尾 len_hinge 个碱基
        - 生成的互补链尾部补 len_hinge 个 'A'
    """
    # 分段
    a_segs = [seq_a[i*seg_len:(i+1)*seg_len] for i in range(3)]

    # hinge 处理函数：取段的前 (seg_len - hinge) 个碱基用于互补
    def hinge_ignored(seg): return seg[:-len_hinge] if len_hinge > 0 else seg
    def hinge_added(seg): return seg + 'A' * len_hinge

    b_segs = [
        seq_b_input[0:seg_len],                        # B[0] - 用户输入（保留原样）
        seq_b_input[seg_len:2*seg_len],                # B[1] - 用户输入（保留原样）
        hinge_added(revc(hinge_ignored(a_segs[0])))    # B[2] - 忽略hinge后互补，尾部补A
    ]

    c_segs = [
        seq_c_input,                                    # C[0] - 用户输入（保留原样）
        hinge_added(revc(hinge_ignored(a_segs[1]))),   # C[1] - 忽略hinge后互补，尾部补A
        hinge_added(revc(hinge_ignored(b_segs[1])))     # C[2] - 忽略hinge后互补，尾部补A
    ]

    d_segs = [
        hinge_added(revc(hinge_ignored(a_segs[2]))),    # D[0] - 忽略hinge后互补，尾部补A
        hinge_added(revc(hinge_ignored(c_segs[0]))),    # D[1] - 忽略hinge后互补，尾部补A
        hinge_added(revc(hinge_ignored(b_segs[0])))     # D[2] - 忽略hinge后互补，尾部补A
    ]

    # 组合成完整链
    seq_b = ''.join(b_segs)
    seq_c = ''.join(c_segs)
    seq_d = ''.join(d_segs)

    # 验证所有链长度一致
    if not (len(seq_a) == len(seq_b) == len(seq_c) == len(seq_d)):
        raise ValueError("补全后链长度不一致")

    # 验证6条边互补性（忽略 hinge 区）
    edges = [
        (hinge_ignored(a_segs[0]), hinge_ignored(b_segs[2]), "A[0]-B[2]"),
        (hinge_ignored(a_segs[1]), hinge_ignored(c_segs[1]), "A[1]-C[1]"),
        (hinge_ignored(a_segs[2]), hinge_ignored(d_segs[0]), "A[2]-D[0]"),
        (hinge_ignored(b_segs[0]), hinge_ignored(d_segs[2]), "B[0]-D[2]"),
        (hinge_ignored(b_segs[1]), hinge_ignored(c_segs[2]), "B[1]-C[2]"),
        (hinge_ignored(c_segs[0]), hinge_ignored(d_segs[1]), "C[0]-D[1]"),
    ]

    print("\n=== 自动补全结果 ===")
    print(f"A链: {seq_a}")
    print(f"B链: {seq_b}")
    print(f"C链: {seq_c}")
    print(f"D链: {seq_d}")
    print()

    all_valid = True
    for seg1, seg2, name in edges:
        is_comp = (seg1 == revc(seg2))
        status = "✓" if is_comp else "✗"
        print(f"  {status} {name} (忽略hinge): {'互补' if is_comp else '不互补'}")
        if not is_comp:
            all_valid = False

    if not all_valid:
        print("\n警告: 部分边不互补，请检查输入序列")
    else:
        print("\n✓ 所有6条边互补性验证通过（hinge区已忽略）")

    return {
        'A': seq_a,
        'B': seq_b,
        'C': seq_c,
        'D': seq_d,
        'len_hinge': len_hinge,
        'seg_len': seg_len,
        'edges_valid': all_valid,
        'segments': {
            'A': a_segs,
            'B': b_segs,
            'C': c_segs,
            'D': d_segs
        }
    }


def generate_standard_file(autofill_data: dict, out_path: str = None):
    """
    将补全后的数据写入标准格式文件

    标准格式:
        len_hinge = N
        A链完整序列
        B链完整序列
        C链完整序列
        D链完整序列
    """
    if out_path is None:
        out_path = 'autofill_input_standard.txt'

    content = f"len_hinge = {autofill_data['len_hinge']}\n"
    content += f"{autofill_data['A']}\n"
    content += f"{autofill_data['B']}\n"
    content += f"{autofill_data['C']}\n"
    content += f"{autofill_data['D']}\n"

    with open(out_path, 'w') as f:
        f.write(content)

    print(f"\n标准格式文件已保存: {out_path}")
    return out_path


def check_input_mismatch(txt_path: str):
    """
    检测输入文件是否为4链标准格式（6段匹配）
    如果检测到4条完整链，提示用户使用标准模式

    返回: (is_standard_format, message)
    """
    with open(txt_path, 'r') as f:
        raw_lines = [line.strip() for line in f]

    if not raw_lines:
        return False, ""

    # 去掉配置行，收集序列
    sequence_lines = raw_lines[1:]

    # 合并多行序列
    sequences_blocks = []
    current_block = []
    for line in sequence_lines:
        if line == '':
            if current_block:
                sequences_blocks.append(''.join(current_block).upper())
                current_block = []
        else:
            current_block.append(line)
    if current_block:
        sequences_blocks.append(''.join(current_block).upper())

    # 如果只有4个非空行，可能是标准格式
    non_empty = [l.upper() for l in sequence_lines if l]
    if len(non_empty) == 4 or len(sequences_blocks) == 4:
        # 检查是否所有4条链等长
        if len(sequences_blocks) == 4:
            lens = [len(s) for s in sequences_blocks]
            if lens[0] == lens[1] == lens[2] == lens[3] and lens[0] % 3 == 0:
                return True, (
                    f"检测到输入文件包含4条等长序列（标准格式），\n"
                    f"每条长度 {lens[0]}，段长 {lens[0]//3}。\n"
                    f"请去掉 -autofill 参数以使用标准模式，\n"
                    f"或使用正确的 -autofill 格式（A完整 + B2段 + C1段）。"
                )

    return False, ""


def process_autofill(txt_path: str, out_path: str = None, verbose: bool = True):
    """
    主入口：处理-autofill模式的输入文件

    1. 检测是否为错误的标准格式
    2. 解析输入
    3. 自动补全序列
    4. 生成标准格式输出文件

    返回: 标准格式文件路径
    """
    # 步骤1: 检测格式
    is_standard, msg = check_input_mismatch(txt_path)
    if is_standard:
        raise ValueError(msg)

    # 步骤2: 解析
    len_hinge, seq_a, seq_b_input, seq_c_input, seg_len = parse_autofill_input(txt_path)

    if verbose:
        print(f"\n[Autofill] 解析输入文件: {txt_path}")
        print(f"[Autofill] len_hinge = {len_hinge}")
        print(f"[Autofill] A链长度: {len(seq_a)} (3段 x {seg_len})")
        print(f"[Autofill] B链输入长度: {len(seq_b_input)} (2段)")
        print(f"[Autofill] C链输入长度: {len(seq_c_input)} (1段)")
        print(f"[Autofill] hinge长度: {len_hinge} (互补匹配时忽略每段末尾{len_hinge}个碱基)")

    # 步骤3: 自动补全
    autofill_data = autofill_sequences(len_hinge, seq_a, seq_b_input, seq_c_input, seg_len)

    # 步骤4: 生成标准文件
    if out_path is None:
        base_name = os.path.splitext(os.path.basename(txt_path))[0]
        out_path = f"{base_name}_autofilled.txt"

    standard_path = generate_standard_file(autofill_data, out_path)

    return standard_path


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python autofill_sequences.py <输入文件.txt> [输出文件.txt]")
        sys.exit(1)

    txt_file = sys.argv[1]
    out_file = sys.argv[2] if len(sys.argv) > 2 else None

    try:
        result = process_autofill(txt_file, out_file)
        print(f"\n成功生成: {result}")
    except ValueError as e:
        print(f"\n错误: {e}")
        sys.exit(1)