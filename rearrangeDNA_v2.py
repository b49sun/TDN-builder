# 输出完整的修改后 rearrangeDNA.py 代码
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Windows GUI 版：互换指定链的两段残基，保持链原始顺序不变
其余链完全原样保留，可多次运行不丢链
"""
import os, sys, tkinter as tk
from tkinter import filedialog, messagebox
from collections import defaultdict

def parse_range(s):
    b, e = map(int, s.split('-'))
    return b, e


def rearrange_chain(infile: str, outfile: str, chain_id: str, cut_position: int):
    """
    对指定链进行残基重排：将 [0:cut_position] 和 [cut_position:] 对换
    
    参数:
        infile: 输入PDB文件路径
        outfile: 输出PDB文件路径  
        chain_id: 要调整的链ID（如 'A', 'B'）
        cut_position: 切割位置，表示修正链[0]在原始链中的位置
                     新顺序 = [cut_position : end] + [0 : cut_position]
                     这样就将修正链恢复为原始链顺序
    """
    # 一次扫描：记录每条链的行
    chain_rows = defaultdict(list)      # 链 -> 原子行
    non_atom = []                       # 非原子行
    
    with open(infile) as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                ch = line[21]
                chain_rows[ch].append(line)
            else:
                non_atom.append(line)

    if chain_id not in chain_rows:
        raise ValueError(f'链 {chain_id} 不存在于PDB文件中')
    
    # 对目标链生成新残基顺序
    res_dict = defaultdict(list)
    for line in chain_rows[chain_id]:
        res_num = int(line[22:26])
        res_dict[res_num].append(line)
    
    all_res = sorted(res_dict.keys())
    n_res = len(all_res)
    
    # cut_position 是修正链[0]在原始链中的位置（从0开始）
    # 重排逻辑：将 [0:cut_position] 移到 [cut_position:] 后面
    # 即新顺序 = [cut_position : n_res] + [0 : cut_position]
    
    if cut_position <= 0 or cut_position >= n_res:
        # 无需重排，直接复制
        with open(outfile, 'w') as fo:
            for ch in sorted(chain_rows.keys()):
                fo.writelines(chain_rows[ch])
            for line in non_atom:
                fo.write(line)
        return True
    
    # 新顺序：后半部分（从cut_position开始）+ 前半部分（cut_position之前）
    new_order = (
        list(range(n_res-cut_position + 1, n_res + 1)) +  # 后半部分（修正链的起始部分）
        list(range(1, n_res-cut_position + 1))           # 前半部分（移到末尾）
    )
    
    # 重写目标链，重新编号残基为1,2,3...
    new_target_lines = []
    res_serial = 1
    for res_num in new_order:
        for line in res_dict[res_num]:
            # 更新残基序号（保持其他信息不变）
            new_line = line[:22] + f"{res_serial:4d}" + line[26:]
            new_target_lines.append(new_line)
        res_serial += 1
    
    chain_rows[chain_id] = new_target_lines

    # 全局原子号重编 & 按原始链顺序写回
    atom_serial = 1
    with open(outfile, 'w') as fo:
        # 按链字母原始顺序输出
        for ch in sorted(chain_rows.keys()):
            for line in chain_rows[ch]:
                new_line = line[:6] + f"{atom_serial:5d}" + line[11:]
                fo.write(new_line)
                atom_serial += 1
        # 非原子行（TER/END/CRYST1...）
        for line in non_atom:
            fo.write(line)
    
    return True


def rearrange_all_chains(infile: str, outfile: str, cut_positions: dict):
    """
    对多条链进行残基重排
    
    参数:
        infile: 输入PDB文件路径
        outfile: 输出PDB文件路径
        cut_positions: {链ID: 切割位置} 字典
                      例如: {'A': 21, 'B': 21, 'C': 21, 'D': 21}
                      切割位置表示修正链[0]在原始链中的索引
    """
    import tempfile
    import shutil
    
    current_file = infile
    temp_files = []
    
    try:
        chains = sorted(cut_positions.keys())
        
        for i, chain_id in enumerate(chains):
            cut_pos = cut_positions[chain_id]
            
            if i == len(chains) - 1:
                # 最后一条链，输出到最终文件
                rearrange_chain(current_file, outfile, chain_id, cut_pos)
            else:
                # 中间步骤，使用临时文件
                fd, temp_path = tempfile.mkstemp(suffix='.pdb')
                os.close(fd)
                rearrange_chain(current_file, temp_path, chain_id, cut_pos)
                temp_files.append(temp_path)
                current_file = temp_path
                
        return True
        
    finally:
        # 清理临时文件
        for temp_path in temp_files:
            if os.path.exists(temp_path):
                os.remove(temp_path)


def main():
    root = tk.Tk()
    root.withdraw()

    infile = filedialog.askopenfilename(
        title='请选择要处理的 PDB 文件',
        filetypes=[('PDB files', '*.pdb'), ('All files', '*.*')])
    if not infile:
        return

    outfile = filedialog.asksaveasfilename(
        title='请指定输出文件',
        defaultextension='.pdb',
        filetypes=[('PDB files', '*.pdb'), ('All files', '*.*')])
    if not outfile:
        return

    target_chain = input('请输入要调整的链ID（如 A、B、C），直接回车跳过：').strip().upper()
    if target_chain:
        seg1_str = input('请输入第一段残基区间（如 1-3）：').strip()
        seg2_str = input('请输入第二段残基区间（如 4-5）：').strip()
        seg1 = parse_range(seg1_str)
        seg2 = parse_range(seg2_str)
        s1_b, s1_e = seg1
        s2_b, s2_e = seg2
        if s1_e >= s2_b:
            messagebox.showerror('错误', '两段区间必须有序且不重叠！')
            return
    else:
        seg1 = seg2 = None

    # 一次扫描：记录每条链的行，以及"非原子行"出现位置
    chain_rows = defaultdict(list)      # 链 -> 原子行
    non_atom = []                       # (行号, 行内容) 用于按位置插回
    line_no = 0
    with open(infile) as f:
        for line in f:
            line_no += 1
            if line.startswith(('ATOM  ', 'HETATM')):
                ch = line[21]
                chain_rows[ch].append(line)
            else:
                non_atom.append((line_no, line))

    # 若跳过，原样复制
    if not target_chain or target_chain not in chain_rows:
        with open(outfile, 'w') as fo:
            for ch in sorted(chain_rows.keys()):
                fo.writelines(chain_rows[ch])
            for _, line in non_atom:
                fo.write(line)
        print('跳过处理，已原样复制。输出文件：', outfile)
        os.system('pause')
        return

    # 对目标链生成新残基顺序
    res_dict = defaultdict(list)
    for line in chain_rows[target_chain]:
        res_num = int(line[22:26])
        res_dict[res_num].append(line)
    all_res = sorted(res_dict.keys())
    new_order = (
        [r for r in all_res if r < seg1[0]] +
        list(range(seg2[0], seg2[1] + 1)) +
        [r for r in all_res if seg1[1] < r < seg2[0]] +
        list(range(seg1[0], seg1[1] + 1)) +
        [r for r in all_res if r > seg2[1]]
    )

    # 重写目标链
    new_target_lines = []
    res_serial = 1
    for res in new_order:
        for line in res_dict[res]:
            new_line = line[:22] + f"{res_serial:4d}" + line[26:]
            new_target_lines.append(new_line)
        res_serial += 1
    chain_rows[target_chain] = new_target_lines

    # 全局原子号重编 & 按原始链顺序写回
    atom_serial = 1
    with open(outfile, 'w') as fo:
        # 按链字母原始顺序输出
        for ch in sorted(chain_rows.keys()):
            for line in chain_rows[ch]:
                new_line = line[:6] + f"{atom_serial:5d}" + line[11:]
                fo.write(new_line)
                atom_serial += 1
        # 非原子行（TER/END/CRYST1...）按原位置写回文件末尾即可
        for _, line in non_atom:
            fo.write(line)

    print('已完成！输出文件：', outfile)
    os.system('pause')

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        messagebox.showerror('运行出错', str(e))
        os.system('pause')
