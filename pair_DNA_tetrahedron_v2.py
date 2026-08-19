# -*- coding: utf-8 -*-
"""
DNA 四面体 V2 版本
输出：修正链（从匹配区开始）和剪切位置（修正链[0]在原始链中的位置）
"""
import numpy as np
from typing import Dict, Tuple, List, Optional

class ValidationError(Exception):
    pass

class TetrahedronCore:
    """
    V2核心类 - 修复版：构建所有6条边
    """
    
    # 所有6条边
    ALL_EDGES = ['EF', 'FG', 'GE', 'EH', 'FH', 'GH']
    # A链的三条边（底部三角形）
    A_EDGES = ['EF', 'FG', 'GE']
    # 顶点对应关系：A段索引 -> (起点, 终点)
    A_VERTICES = [
        ('E', 'F'),  # A[0]: EF边
        ('F', 'G'),  # A[1]: FG边
        ('G', 'E')   # A[2]: GE边
    ]
    
    def __init__(self, sequences: Dict[str, str], ignore_last_bases: int = 0, verbose: bool = True):
        if len(sequences) != 4:
            raise ValidationError('必须提供4条链')
        
        # 找到A链
        if 'A' in sequences:
            self.seq_a_original = sequences['A'].upper()
            candidate_seqs = {k: v.upper() for k, v in sequences.items() if k != 'A'}
        else:
            keys = list(sequences.keys())
            self.seq_a_original = sequences[keys[0]].upper()
            candidate_seqs = {k: sequences[k].upper() for k in keys[1:]}
        
        self.L = len(self.seq_a_original)
        self.seg_len = self.L // 3
        self.len_hinge = ignore_last_bases if ignore_last_bases > 0 else 1
        self.ignore_last_bases = ignore_last_bases
        self.verbose = verbose
        
        if self.L % 3 != 0:
            raise ValidationError('链长必须是3的倍数')
        if any(len(s) != self.L for s in candidate_seqs.values()):
            raise ValidationError('所有链必须等长')
        
        self.candidate_seqs_original = candidate_seqs
        self.comp = str.maketrans('ATGC', 'TACG')
        
        # 为候选链创建环形扩展序列
        self.extended_seqs = {}
        for ch, seq in candidate_seqs.items():
            self.extended_seqs[ch] = seq + seq[:self.seg_len]
    
    def _revc(self, s: str) -> str:
        return s.translate(self.comp)[::-1]
    
    def _get_edge_name(self, v1: str, v2: str) -> str:
        """无向边名称（排序）"""
        return ''.join(sorted([v1, v2]))
    
    def find_complementary_segment(self, target_seq: str, search_ch: str, head_offset: int = 0) -> Optional[Tuple[int, int]]:
        """
        寻找互补段
        head_offset: 头部hinge碱基数，用于调整匹配起始位置
                    hinge=1时为0，hinge=3时为1
        """
        extended = self.extended_seqs[search_ch]
    
        if target_seq not in extended:
            return None
    
        start = extended.index(target_seq)
        # 调整起始位置：往前推head_offset个碱基，包含头部hinge
        adjusted_start = (start - head_offset) % self.L
        orig_start = adjusted_start % self.L
    
        return (orig_start, len(target_seq))
    
    def try_partition(self, cut1: int) -> Optional[Dict]:
        """尝试A链在cut1处开口 - 强制使用 head_tail 模式"""
        seq_a = self.seq_a_original
        L = self.L
        seg_len = self.seg_len
        
        # 强制使用 head_tail 模式
        if self.ignore_last_bases == 1:
            head, tail = 0, 1
        elif self.ignore_last_bases == 2:
            head, tail = 1, 1  # 均分
        elif self.ignore_last_bases == 3:
            head, tail = 1, 2
        elif self.ignore_last_bases == 4:
            head, tail = 2, 2
        else:
            head, tail = 0, self.ignore_last_bases
        
        # A链剪切位置往后推head个碱基
        adjusted_cut1 = (cut1 + head) % L
        
        # A链三段（从adjusted_cut1开始的循环序列）
        new_seq = seq_a[adjusted_cut1:] + seq_a[:adjusted_cut1]
        seg_a0_full = new_seq[0:seg_len]
        seg_a1_full = new_seq[seg_len:2*seg_len]
        seg_a2_full = new_seq[2*seg_len:]
        
        # 用于匹配的序列（去掉头部和尾部hinge区）
        if head > 0 or tail > 0:
            seg_a0 = seg_a0_full[head:-tail] if len(seg_a0_full) > (head + tail) else ""
            seg_a1 = seg_a1_full[head:-tail] if len(seg_a1_full) > (head + tail) else ""
            seg_a2 = seg_a2_full[head:-tail] if len(seg_a2_full) > (head + tail) else ""
        else:
            seg_a0, seg_a1, seg_a2 = seg_a0_full, seg_a1_full, seg_a2_full
        
        # 为每段寻找互补链
        available = list(self.candidate_seqs_original.keys())
        assignment = {}
        dynamic_names = {}
        
        for a_idx in range(3):
            target_rev = self._revc([seg_a0, seg_a1, seg_a2][a_idx])
            found = False
            
            for ch in available[:]:
                result = self.find_complementary_segment(target_rev, ch, head_offset=head)
                if result is not None:
                    match_start, match_len = result
                    assignment[a_idx] = (ch, match_start, match_len)
                    dynamic_names[a_idx] = ['B', 'C', 'D'][a_idx]
                    available.remove(ch)
                    found = True
                    break
            
            if not found:
                return None
        
        # ==================== 新增验证：B/C/D链之间的三条边 ====================
        
        # 构建B/C/D链的段信息
        bcd_info = {}  # {动态链名: 详细信息}
        
        for a_idx, (orig_ch, match_start, match_len) in assignment.items():
            seq_original = self.candidate_seqs_original[orig_ch]
            # 修正链：从match_start开始
            corrected_seq = seq_original[match_start:] + seq_original[:match_start]
            dyn_name = dynamic_names[a_idx]
            
            # 修正链的三段（每段包含hinge）
            seg0 = corrected_seq[0:seg_len]
            seg1 = corrected_seq[seg_len:2*seg_len]
            seg2 = corrected_seq[2*seg_len:]
            
            # 用于匹配的三段（去掉hinge）
            if head > 0 or tail > 0:
                match_seg0 = seg0[head:-tail] if len(seg0) > (head + tail) else ""
                match_seg1 = seg1[head:-tail] if len(seg1) > (head + tail) else ""
                match_seg2 = seg2[head:-tail] if len(seg2) > (head + tail) else ""
            else:
                match_seg0, match_seg1, match_seg2 = seg0, seg1, seg2
            
            # 确定与A匹配的段索引
            # match_start是原始链中的位置，需要映射到修正链的段索引
            # 修正链从match_start开始，所以原始链的match_start位置现在是修正链的0位置
            # 因此，与A匹配的是修正链的段0
            matched_idx_in_corrected = 0  # 修正链的段0与A匹配
            
            # 计算另外两段在原始链中的位置
            other_segs = []
            for idx in range(3):
                if idx == matched_idx_in_corrected:
                    continue
                start_in_corrected = idx * seg_len
                end_in_corrected = (idx + 1) * seg_len
                # 映射回原始链
                start_in_original = (match_start + start_in_corrected) % L
                end_in_original = (match_start + end_in_corrected) % L
                other_segs.append({
                    'corrected_idx': idx,
                    'start_orig': start_in_original,
                    'end_orig': end_in_original,
                    'match_seq': [match_seg0, match_seg1, match_seg2][idx]
                })
            
            bcd_info[dyn_name] = {
                'orig_ch': orig_ch,
                'match_start': match_start,
                'corrected_seq': corrected_seq,
                'full_segs': [seg0, seg1, seg2],
                'match_segs': [match_seg0, match_seg1, match_seg2],
                'matched_idx': matched_idx_in_corrected,
                'other_segs': other_segs
            }
        
        # 验证B/C/D之间的三条边
        bcd_names = ['B', 'C', 'D']
        
        for i, name in enumerate(bcd_names):
            current = bcd_info[name]
            other_names = [n for n in bcd_names if n != name]
            
            for seg_info in current['other_segs']:
                idx = seg_info['corrected_idx']
                match_seq = seg_info['match_seq']
                target_rev = self._revc(match_seq)
                
                # 检查是否能与另外两条链的某段互补
                found_match = False
                for other_name in other_names:
                    other_info = bcd_info[other_name]
                    for other_idx, other_match_seq in enumerate(other_info['match_segs']):
                        if other_idx == other_info['matched_idx']:
                            continue  # 跳过与A匹配的段
                        if other_match_seq == target_rev:
                            found_match = True
                            break
                    if found_match:
                        break
                
                if not found_match:
                    #print(f"    -> 未找到匹配!")
                    return None
        
        # 构建修正链和剪切位置
        corrected_seqs = {'A': new_seq}
        cut_positions = {'A': adjusted_cut1}
        
        for a_idx, (orig_ch, match_start, match_len) in assignment.items():
            seq_original = self.candidate_seqs_original[orig_ch]
            corrected_seq = seq_original[match_start:] + seq_original[:match_start]
            dyn_name = dynamic_names[a_idx]
            corrected_seqs[dyn_name] = corrected_seq
            cut_positions[dyn_name] = match_start
        
        cut2 = (adjusted_cut1 + seg_len) % L
        return {
            'cuts': (adjusted_cut1, cut2),
            'assignment': assignment,
            'dynamic_names': dynamic_names,
            'corrected_sequences': corrected_seqs,
            'cut_positions': cut_positions
        }
    
    def find_all_valid_partitions(self) -> List[Dict]:
        """遍历cut1从0到seg_len-1"""
        valid = []
        seg_len = self.seg_len
        
        for cut1 in range(0, seg_len):
            r = self.try_partition(cut1)
            if r:
                valid.append(r)
        
        return valid
    
    def validate_edges(self) -> Dict:
        """验证边配对 - 修复版：构建所有6条边"""
        partitions = self.find_all_valid_partitions()
        
        if not partitions:
            raise ValidationError('未找到有效的四面体配置，请检查序列互补性')
        
        best = partitions[0]
        self.best_partition = best
        
        # 打印修正链和剪切位置
        if self.verbose:
            for ch in ['A', 'B', 'C', 'D']:
                seq = best['corrected_sequences'][ch]
                cut_pos = best['cut_positions'][ch]
        
        # 构建边配对 - 修复：构建所有6条边
        assignment = best['assignment']
        dynamic_names = best['dynamic_names']
        
        # 初始化所有6条边
        edge_data = {name: [] for name in self.ALL_EDGES}
        a_to_other = {}  # {A段索引: (动态链名, 段索引)}
        
        # 步骤1: 构建A链的三条边（EF, FG, GE）
        for a_idx, (orig_ch, match_start, match_len) in assignment.items():
            edge_name = self.A_EDGES[a_idx]
            dyn_name = dynamic_names[a_idx]
            seg_idx = match_start // self.seg_len
            
            edge_data[edge_name].append(('A', a_idx))
            edge_data[edge_name].append((dyn_name, seg_idx))
            a_to_other[a_idx] = (dyn_name, seg_idx)
        
        # 步骤2: 为B/C/D构建连接到H的边（EH, FH, GH）
        for a_idx, (dyn_name, matched_seg_idx) in a_to_other.items():
            a_start, a_end = self.A_VERTICES[a_idx]
            # 反向：A链从start到end，互补链从end到start
            match_start_vertex, match_end_vertex = a_end, a_start
            
            # 计算另外两个段的位置（循环）
            seg_next = (matched_seg_idx + 1) % 3
            seg_prev = (matched_seg_idx - 1) % 3
            
            # 连接到H的边
            edge_next = self._get_edge_name(match_end_vertex, 'H')
            edge_prev = self._get_edge_name('H', match_start_vertex)
            
            edge_data[edge_next].append((dyn_name, seg_next))
            edge_data[edge_prev].append((dyn_name, seg_prev))
        
        # 步骤3: 验证所有边都有2段
        for edge_name, segments in edge_data.items():
            if len(segments) != 2:
                raise ValidationError(f'边{edge_name}应有2段，实际{len(segments)}段')
        
        # 步骤4: 构建输出格式
        edge_assignment = {}
        for edge_name, segments in edge_data.items():
            edge_assignment[edge_name] = {
                'chains': (segments[0][0], segments[1][0]),
                'segments': (segments[0][1], segments[1][1])
            }
        
        return edge_assignment
    
    def get_corrected_data(self) -> Tuple[Dict, Dict]:
        """获取修正链和剪切位置"""
        if hasattr(self, 'best_partition'):
            p = self.best_partition
            return p['corrected_sequences'], p['cut_positions']
        raise ValidationError('请先调用validate_edges()')


def find_valid_configuration(sequences: Dict[str, str], len_hinge: int = 1, verbose: bool = True):
    """
    模块级便捷函数
    返回: (partition_info, corrected_sequences, cut_positions)
    """
    core = TetrahedronCore(sequences, ignore_last_bases=len_hinge, verbose=verbose)
    core.validate_edges()
    corrected_seqs, cut_positions = core.get_corrected_data()
    return core.best_partition, corrected_seqs, cut_positions