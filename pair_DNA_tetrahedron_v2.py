# -*- coding: utf-8 -*-
"""
DNA 四面体 V2 版本：最小改动回溯匹配
输出：修正链（从匹配区开始）和剪切位置（修正链[0]在原始链中的位置）
"""
import numpy as np
from typing import Dict, Tuple, List, Optional

class ValidationError(Exception):
    pass

class TetrahedronCore:
    """
    V2核心类：保留V2拓扑规则，增加匹配回溯搜索
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
    
    def find_complementary_segments(self, target_seq: str, search_ch: str, head_offset: int = 0) -> List[Tuple[int, int]]:
        """Return every circular match position for a target sequence."""
        extended = self.extended_seqs[search_ch]
        matches = []
        seen = set()
        # A valid circular start is restricted to the original strand (0..L-1).
        for start in range(self.L):
            if extended[start:start + len(target_seq)] == target_seq:
                adjusted_start = (start - head_offset) % self.L
                item = (adjusted_start, len(target_seq))
                if item not in seen:
                    matches.append(item)
                    seen.add(item)
        return matches

    def find_complementary_segment(self, target_seq: str, search_ch: str, head_offset: int = 0) -> Optional[Tuple[int, int]]:
        """Backward-compatible helper: return the first match, if any."""
        matches = self.find_complementary_segments(target_seq, search_ch, head_offset)
        return matches[0] if matches else None

    def _hinge_offsets(self) -> Tuple[int, int]:
        if self.ignore_last_bases == 1:
            return 0, 1
        if self.ignore_last_bases == 2:
            return 1, 1
        if self.ignore_last_bases == 3:
            return 1, 2
        if self.ignore_last_bases == 4:
            return 2, 2
        return 0, self.ignore_last_bases

    def _candidate_passes_edge_checks(self, assignment: Dict[int, Tuple[str, int, int]],
                                      dynamic_names: Dict[int, str],
                                      head: int, tail: int,
                                      a_sequence: str) -> bool:
        """Validate the six *designated* tetrahedral edges for this candidate.

        The corrected chain for each B/C/D candidate starts at its A-matching
        segment. Therefore segment 0 of each corrected B/C/D chain must pair
        with the corresponding A segment. The remaining four edge pairings
        are fixed by the tetrahedral strand connectivity:

            EF: A[0] -- B[0]       EH: B[1] -- D[2]
            FG: A[1] -- C[0]       FH: B[2] -- C[1]
            GE: A[2] -- D[0]       GH: C[2] -- D[1]

        This replaces the former permissive "matches any other segment" test.
        It deliberately checks the exact partner prescribed by the topology.
        """
        seg_len = self.seg_len

        # Reconstruct the exact corrected sequences associated with this
        # candidate assignment. Do not infer edge partners from arbitrary
        # sequence matches elsewhere in the other chains.
        corrected = {'A': a_sequence}
        for a_idx, (orig_ch, match_start, _match_len) in assignment.items():
            original = self.candidate_seqs_original[orig_ch]
            corrected[dynamic_names[a_idx]] = (
                original[match_start:] + original[:match_start]
            )

        if set(corrected) != {'A', 'B', 'C', 'D'}:
            return False
        if any(len(seq) != self.L for seq in corrected.values()):
            return False

        # Split each corrected strand into its three consecutive segments,
        # then remove the hinge nucleotides using the same convention as the
        # A-chain matching step.
        segments = {}
        for chain_name, sequence in corrected.items():
            full_segments = [
                sequence[:seg_len],
                sequence[seg_len:2 * seg_len],
                sequence[2 * seg_len:3 * seg_len],
            ]
            if head > 0 or tail > 0:
                match_segments = [
                    segment[head:-tail]
                    if len(segment) > head + tail else ''
                    for segment in full_segments
                ]
            else:
                match_segments = full_segments

            if any(not segment for segment in match_segments):
                return False
            segments[chain_name] = match_segments

        # Exact edge-to-segment topology. Each tuple identifies
        # (chain, segment_index) for the two strands forming that edge.
        expected_edges = {
            'EF': (('A', 0), ('B', 0)),
            'FG': (('A', 1), ('C', 0)),
            'GE': (('A', 2), ('D', 0)),
            'EH': (('B', 1), ('D', 2)),
            'FH': (('B', 2), ('C', 1)),
            'GH': (('C', 2), ('D', 1)),
        }

        # Every edge must contain exactly its designated pair, and the two
        # segments must be exact reverse complements of one another.
        for edge_name, ((chain_1, index_1), (chain_2, index_2)) in expected_edges.items():
            seq_1 = segments[chain_1][index_1]
            seq_2 = segments[chain_2][index_2]
            if len(seq_1) != len(seq_2):
                return False
            if self._revc(seq_1) != seq_2:
                return False

        # Retain an explicit six-edge occupancy sanity check. This guards
        # against accidental omissions if the topology table is edited later.
        edge_data = {edge_name: [] for edge_name in self.ALL_EDGES}
        for edge_name, pair in expected_edges.items():
            edge_data[edge_name].extend(pair)
        return all(len(pair) == 2 for pair in edge_data.values())

    def try_partition(self, cut1: int) -> Optional[Dict]:
        """Try candidate matches with backtracking; return the first fully valid one."""
        seq_a, L, seg_len = self.seq_a_original, self.L, self.seg_len
        head, tail = self._hinge_offsets()
        adjusted_cut1 = (cut1 + head) % L
        new_seq = seq_a[adjusted_cut1:] + seq_a[:adjusted_cut1]
        full_a = [new_seq[:seg_len], new_seq[seg_len:2 * seg_len], new_seq[2 * seg_len:]]
        if head > 0 or tail > 0:
            match_a = [s[head:-tail] if len(s) > head + tail else '' for s in full_a]
        else:
            match_a = full_a
        if any(not s for s in match_a):
            return None

        # Collect all candidate matches for each A segment before assigning chains.
        options = []
        for a_idx in range(3):
            target_rev = self._revc(match_a[a_idx])
            candidates = []
            for ch in self.candidate_seqs_original:
                for match_start, match_len in self.find_complementary_segments(
                        target_rev, ch, head_offset=head):
                    candidates.append((ch, match_start, match_len))
            if not candidates:
                return None
            options.append(candidates)

        dynamic_names = {0: 'B', 1: 'C', 2: 'D'}
        assignment = {}

        def search(a_idx: int, used_chains: set):
            if a_idx == 3:
                if not self._candidate_passes_edge_checks(assignment, dynamic_names, head, tail, new_seq):
                    return None

                corrected_seqs = {'A': new_seq}
                cut_positions = {'A': adjusted_cut1}
                for idx, (orig_ch, match_start, _) in assignment.items():
                    original = self.candidate_seqs_original[orig_ch]
                    corrected_seqs[dynamic_names[idx]] = original[match_start:] + original[:match_start]
                    cut_positions[dynamic_names[idx]] = match_start

                return {
                    'cuts': (adjusted_cut1, (adjusted_cut1 + seg_len) % L),
                    'assignment': dict(assignment),
                    'dynamic_names': dict(dynamic_names),
                    'corrected_sequences': corrected_seqs,
                    'cut_positions': cut_positions
                }

            for ch, match_start, match_len in options[a_idx]:
                if ch in used_chains:
                    continue
                assignment[a_idx] = (ch, match_start, match_len)
                result = search(a_idx + 1, used_chains | {ch})
                if result is not None:
                    return result
                del assignment[a_idx]
            return None

        return search(0, set())

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
