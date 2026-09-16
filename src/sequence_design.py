"""
小核酸药物序列设计模块
包含siRNA、miRNA、ASO等小核酸的设计和优化功能
包含四套经典siRNA评分规则：Reynolds、Ui-Tei、Amarzguioui、热力学不对称性
"""

import re
import math
from typing import List, Tuple, Optional, Dict

class NucleicAcidDesigner:
    """
    小核酸序列设计器
    """

    def __init__(self):
        self.complement_map = {
            'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C',
            'a': 'u', 't': 'a', 'c': 'g', 'g': 'c',
            'U': 'A', 'u': 'a'
        }

        self.nn_params = {
            'AA': -1.0, 'AU': -0.9, 'UA': -1.1, 'UU': -0.9,
            'AG': -1.3, 'GA': -1.4, 'CG': -2.1, 'GC': -2.7,
            'GG': -1.8, 'CA': -1.4, 'AC': -1.5, 'CC': -1.8,
            'UC': -1.3, 'CU': -1.3, 'GU': -1.4, 'UG': -1.5,
        }

        self.immune_stimulatory_sequences = [
            ('UGUGU', 'TLR7/8', '免疫激活'),
            ('GUCCUUCAA', 'TLR7', '免疫激活'),
        ]

    def complement(self, sequence: str) -> str:
        sequence = sequence.replace('T', 'U').replace('t', 'u')
        return ''.join([self.complement_map.get(base, base) for base in sequence])

    def reverse_complement(self, sequence: str) -> str:
        return self.complement(sequence)[::-1]

    def calculate_gc_content(self, sequence: str) -> float:
        gc_count = sequence.count('G') + sequence.count('C') + sequence.count('g') + sequence.count('c')
        return gc_count / len(sequence) * 100

    def is_valid_sequence(self, sequence: str) -> bool:
        valid_bases = {'A', 'T', 'C', 'G', 'U', 'a', 't', 'c', 'g', 'u', 'N', 'n'}
        return all(base in valid_bases for base in sequence)

    def calculate_delta_g(self, sequence: str) -> float:
        if len(sequence) < 2:
            return 0.0
        total_dg = 0.0
        seq_upper = sequence.upper().replace('T', 'U')
        for i in range(len(seq_upper) - 1):
            dinuc = seq_upper[i:i+2]
            total_dg += self.nn_params.get(dinuc, -1.0)
        return total_dg

    def get_5prime_stability(self, sequence: str, length: int = 4) -> float:
        return abs(self.calculate_delta_g(sequence[:length]))

    def design_sirna(self, target_mrna: str, position: int = 100, length: int = 21) -> Tuple[str, str]:
        if position + length > len(target_mrna):
            raise ValueError("目标位置超出mRNA序列范围")
        sense_strand = target_mrna[position-1 : position-1+length]
        antisense_strand = self.reverse_complement(sense_strand)
        return sense_strand, antisense_strand

    def check_immune_stimulatory(self, sense: str) -> Tuple[bool, List[str]]:
        """
        检查免疫刺激序列 (HF4)
        """
        sense = sense.upper().replace('T', 'U')
        warnings = []
        
        for seq, target, risk in self.immune_stimulatory_sequences:
            if seq in sense:
                warnings.append(f"含{seq}({target}): {risk}")
        
        if re.search(r'UUUUUU', sense):
            warnings.append("含≥6连续U: TLR信号传导风险")
        
        return len(warnings) == 0, warnings

    def blast_alignment(self, query: str, subject: str) -> Tuple[float, str]:
        """
        BLAST比对算法 - 种子区匹配 + 全局比对
        考虑T/U转换
        :param query: 查询序列
        :param subject: 目标序列
        :return: (相似度百分比, 比对详情)
        """
        query = query.upper().replace('T', 'U')
        subject = subject.upper().replace('T', 'U')
        
        query_len = len(query)
        subject_len = len(subject)
        
        if query_len == 0 or subject_len == 0:
            return 0.0, "序列为空"
        
        seed_size = 7
        seed_matches = []
        
        for i in range(query_len - seed_size + 1):
            seed = query[i:i+seed_size]
            for j in range(subject_len - seed_size + 1):
                if seed == subject[j:j+seed_size]:
                    seed_matches.append((i, j))
        
        if not seed_matches:
            return 0.0, "无种子区匹配"
        
        best_score = 0
        best_align = ""
        
        for seed_i, seed_j in seed_matches:
            start_i = max(0, seed_i - 5)
            start_j = max(0, seed_j - 5)
            end_i = min(query_len, seed_i + seed_size + 5)
            end_j = min(subject_len, seed_j + seed_size + 5)
            
            query_segment = query[start_i:end_i]
            subject_segment = subject[start_j:end_j]
            
            matches = sum(1 for a, b in zip(query_segment, subject_segment) if a == b)
            length = max(len(query_segment), len(subject_segment))
            score = matches / length * 100
            
            if score > best_score:
                best_score = score
                best_align = f"位置{start_j+1}-{end_j}: {matches}/{length}"
        
        global_matches = sum(1 for a, b in zip(query, subject[:query_len]) if a == b)
        global_score = global_matches / query_len * 100
        
        final_score = max(best_score, global_score)
        
        return final_score, best_align
    
    def reverse_complement(self, sequence: str) -> str:
        """获取序列的反向互补链"""
        complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join([complement.get(base, base) for base in sequence[::-1]])

    def check_off_target(self, sense: str, mrna_sequence: str) -> Dict:
        """
        检查脱靶效应 - 双链BLAST比对（正义链+反义链）
        :param sense: siRNA正义链
        :param mrna_sequence: 完整mRNA序列
        :return: 脱靶检查结果
        """
        sense = sense.upper().replace('T', 'U')
        antisense = self.reverse_complement(sense)
        mrna = mrna_sequence.upper().replace('T', 'U')

        results = {
            'has_off_target': False,
            'matches': [],
            'max_identity': 0.0,
            'threshold': 80.0
        }

        window_size = len(sense)
        mrna_len = len(mrna)

        if mrna_len < window_size:
            return results

        for i in range(mrna_len - window_size + 1):
            target_segment = mrna[i:i+window_size]

            sense_matches = sum(1 for a, b in zip(sense, target_segment) if a == b)
            sense_identity = (sense_matches / window_size) * 100

            antisense_matches = sum(1 for a, b in zip(antisense, target_segment) if a == b)
            antisense_identity = (antisense_matches / window_size) * 100

            max_identity = max(sense_identity, antisense_identity)
            
            if max_identity >= 80.0 and max_identity < 100.0:
                strand = "正义链" if sense_identity >= antisense_identity else "反义链"
                matches = sense_matches if sense_identity >= antisense_identity else antisense_matches
                
                results['has_off_target'] = True
                results['matches'].append({
                    'position': i + 1,
                    'identity': max_identity,
                    'strand': strand,
                    'alignment': f"{strand} 位置{i+1}-{i+window_size}: {matches}/{window_size}"
                })
                if max_identity > results['max_identity']:
                    results['max_identity'] = max_identity

        return results
    
    def hard_filter(self, sense: str) -> Tuple[bool, Dict]:
        """
        硬性过滤 (4条)
        HF1: GC含量 25-65%
        HF2: 无≥5连续相同碱基
        HF3: Reynolds评分≥2
        HF4: 无免疫刺激序列
        """
        sense = sense.upper().replace('T', 'U')
        gc = self.calculate_gc_content(sense)
        reynolds = self.reynolds_score(sense)
        immune_safe, immune_warnings = self.check_immune_stimulatory(sense)

        hf1_pass = 25.0 <= gc <= 65.0
        hf2_pass = not re.search(r'(.)\1{4,}', sense)
        hf3_pass = reynolds >= 2
        hf4_pass = immune_safe

        all_pass = hf1_pass and hf2_pass and hf3_pass and hf4_pass

        filter_results = {
            'hf1_gc': {'pass': hf1_pass, 'value': gc, 'standard': '25-65%'},
            'hf2_poly': {'pass': hf2_pass, 'value': '无≥5连续', 'standard': '≤5连续'},
            'hf3_reynolds': {'pass': hf3_pass, 'value': reynolds, 'standard': '≥2'},
            'hf4_immune': {'pass': hf4_pass, 'value': immune_warnings, 'standard': '无免疫刺激'},
            'all_pass': all_pass
        }

        return all_pass, filter_results

    def reynolds_score(self, sense: str) -> int:
        """
        Reynolds评分 (2004)
        满分: 8分
        R1: pos19=G (+1)
        R2: pos3=A (+1)
        R3: pos10=U (+1)
        R4: pos1=A/U (+1)
        R5: GC 30-55% (+1)
        R6: 无≥4连续相同碱基 (+1)
        """
        score = 0
        sense = sense.upper().replace('T', 'U')

        if sense[18] == 'G':
            score += 1
        if sense[2] == 'A':
            score += 1
        if sense[9] == 'U':
            score += 1
        if sense[0] in ['A', 'U']:
            score += 1

        gc = self.calculate_gc_content(sense)
        if 30.0 <= gc <= 55.0:
            score += 1

        if not re.search(r'(.)\1{3,}', sense):
            score += 1

        return score

    def amarzguioui_score(self, sense: str) -> int:
        """
        Amarzguioui评分规则 (满分5分)
        AM1: pos1=A/C (+1)
        AM2: pos6=A (+1)
        AM3: pos19≠G (+1)
        AM4: GC 30-55% (+1)
        AM5: 无≥3连续U (+1)
        """
        score = 0
        sense = sense.upper().replace('T', 'U')

        if sense[0] in ['A', 'C']:
            score += 1
        if sense[5] == 'A':
            score += 1
        if sense[18] != 'G':
            score += 1

        gc = self.calculate_gc_content(sense)
        if 30.0 <= gc <= 55.0:
            score += 1

        if not re.search(r'UUU', sense):
            score += 1

        return score

    def uitei_check(self, sense: str) -> Tuple[bool, str]:
        """
        Ui-Tei规则 (Ui-Tei et al., 2004) - 需全部通过
        
        UT1: 正义链5'端第1位 = A/U (低稳定性)
        UT2: 反义链5'端第1位 = G/C (高稳定性)
        UT3: GC含量 30-55%
        UT4: 无≥4个连续相同碱基
        UT5: 正义链5'端前4bp的ΔG > -5.0 kcal/mol (适度不稳定)
        
        使用最近邻热力学参数计算ΔG
        """
        s = sense.upper().replace('T', 'U')
        a = self.reverse_complement(s)

        # UT1: 正义链1位 = A或U
        ut1 = s[0] in ('A', 'U')
        
        # UT2: 反义链1位 = G或C
        ut2 = a[0] in ('G', 'C')
        
        # UT3: GC含量 30-55%
        gc = self.calculate_gc_content(s)
        ut3 = 30.0 <= gc <= 55.0
        
        # UT4: 无≥4连续相同碱基
        has_poly4 = re.search(r'(.)\1{3,}', s)
        ut4 = not has_poly4
        
        # UT5: 正义链5'端ΔG > -5.0 kcal/mol
        s5 = s[:4]
        dg = sum(self.nn_params.get(s5[i:i+2], -1.0) for i in range(len(s5) - 1))
        ut5 = dg > -5.0

        all_pass = ut1 and ut2 and ut3 and ut4 and ut5

        if all_pass:
            return True, "完全符合Ui-Tei规则"
        else:
            reasons = []
            if not ut1:
                reasons.append(f"正义链1位为{s[0]}，应为A/U")
            if not ut2:
                reasons.append(f"反义链1位为{a[0]}，应为G/C")
            if not ut3:
                reasons.append(f"GC含量{gc:.1f}%，不在30-55%")
            if not ut4:
                reasons.append(f"含≥4连续碱基: {has_poly4.group()}")
            if not ut5:
                reasons.append(f"正义链5'端ΔG={dg:.1f}≤-5.0")
            return False, "; ".join(reasons)

    def thermodynamic_asymmetry(self, sense: str) -> Tuple[bool, float]:
        """
        热力学不对称性 (Schwarz规则)
        条件: 反义链5'端应更稳定(ΔG更负)，有利于RISC选择反义链作为引导链
        计算: 比较双链两端4bp的稳定性
        ΔG越负表示越稳定
        """
        sense = sense.upper().replace('T', 'U')
        antisense = self.reverse_complement(sense)

        def calculate_dimer_stability(seq):
            """计算序列的二核苷酸堆叠能总和"""
            if len(seq) < 2:
                return 0.0
            total_dg = 0.0
            for i in range(len(seq) - 1):
                dinuc = seq[i:i+2]
                total_dg += self.nn_params.get(dinuc, -1.0)
            return total_dg

        sense_5prime_dg = calculate_dimer_stability(sense[:4])
        antisense_5prime_dg = calculate_dimer_stability(antisense[:4])

        return antisense_5prime_dg < sense_5prime_dg, antisense_5prime_dg - sense_5prime_dg

    def calculate_seed_complexity(self, sense: str) -> float:
        """
        种子区复杂度 (反义链2-8位)
        计算归一化香农熵: H = -Σ(p_i × log₂(p_i)) / log₂(4)
        """
        antisense = self.reverse_complement(sense.upper().replace('T', 'U'))
        seed = antisense[1:8]
        
        base_count = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
        for base in seed:
            if base in base_count:
                base_count[base] += 1

        total = len(seed)
        if total == 0:
            return 0.0

        entropy = 0.0
        for count in base_count.values():
            if count > 0:
                p = count / total
                entropy -= p * math.log2(p)

        max_entropy = math.log2(4)
        complexity = entropy / max_entropy if max_entropy > 0 else 0.0

        return complexity

    def get_target_region_bonus(self, position: int, mrna_length: int) -> Tuple[int, str]:
        """
        靶区域位置加分
        3'-UTR: +10
        CDS: +5
        5'-UTR: +0
        起始密码子附近: -5
        """
        cds_start = int(mrna_length * 0.25)
        cds_end = int(mrna_length * 0.75)
        
        start_codon_region_start = max(0, cds_start - 50)
        start_codon_region_end = cds_start + 50
        
        if start_codon_region_start <= position <= start_codon_region_end:
            return -5, "起始密码子附近"
        elif position > cds_end:
            return 10, "3'-UTR"
        elif position >= cds_start:
            return 5, "CDS"
        else:
            return 0, "5'-UTR"
    
    def get_target_region_bonus_by_name(self, region_name: str) -> Tuple[int, str]:
        """
        根据靶区域名称直接返回加分
        """
        bonus_map = {
            "3'-UTR": (10, "3'-UTR"),
            "5'-UTR": (0, "5'-UTR"),
            "CDS": (5, "CDS"),
            "起始密码子附近": (-5, "起始密码子附近"),
            "3'UTR": (10, "3'-UTR"),
            "5'UTR": (0, "5'-UTR"),
        }
        return bonus_map.get(region_name, (5, "CDS"))

    def comprehensive_sirna_evaluation(self, sense: str, position: int = 0, mrna_length: int = 0, target_region_name: str = None) -> Dict:
        """
        综合siRNA评估

        参数:
            sense: 正义链序列
            position: mRNA位置（可选）
            mrna_length: mRNA长度（可选）
            target_region_name: 靶区域名称（可选，如"3'-UTR"、"CDS"等）

        如果提供target_region_name，则使用该名称直接获取靶区域加分
        否则使用position和mrna_length动态计算
        """
        sense = sense.upper().replace('T', 'U')
        antisense = self.reverse_complement(sense)
        gc = self.calculate_gc_content(sense)

        hard_pass, hard_results = self.hard_filter(sense)

        reynolds = self.reynolds_score(sense)
        amarzguioui = self.amarzguioui_score(sense)
        seed_complexity = self.calculate_seed_complexity(sense)
        thermo_pass, thermo_diff = self.thermodynamic_asymmetry(sense)
        uitei_pass, uitei_reason = self.uitei_check(sense)

        optimal_gc = 30.0 <= gc <= 55.0

        if target_region_name is not None:
            target_bonus, target_region = self.get_target_region_bonus_by_name(target_region_name)
        else:
            target_bonus, target_region = self.get_target_region_bonus(position, mrna_length)

        total_score = 0.0

        if hard_pass:
            total_score += 30

        total_score += reynolds * 4
        total_score += amarzguioui * 3
        total_score += seed_complexity * 8

        if thermo_pass:
            total_score += 3
        if uitei_pass:
            total_score += 3
        if optimal_gc:
            total_score += 5

        total_score += target_bonus

        evaluation = {
            'sense': sense,
            'antisense': antisense,
            'gc_content': gc,
            'position': position,
            'target_region': target_region,
            'hard_filter': hard_results,
            'passes_hard_filters': hard_pass,
            'reynolds_score': reynolds,
            'amarzguioui_score': amarzguioui,
            'seed_complexity': seed_complexity,
            'thermo_asymmetry': thermo_pass,
            'thermo_diff': thermo_diff,
            'uitei_passed': uitei_pass,
            'uitei_reason': uitei_reason,
            'optimal_gc': optimal_gc,
            'target_bonus': target_bonus,
            'final_score': total_score,
            'max_score': 96,
            'recommendation': ''
        }

        if hard_pass and total_score >= 70:
            evaluation['recommendation'] = '推荐'
        elif hard_pass and total_score >= 50:
            evaluation['recommendation'] = '良好'
        elif hard_pass:
            evaluation['recommendation'] = '一般'
        else:
            evaluation['recommendation'] = '不推荐'

        return evaluation

    def find_sirna_candidates(self, target_mrna: str, start_pos: int, end_pos: int, length: int = 21, max_candidates: int = 20) -> List[Dict]:
        """
        在指定区域内查找最佳siRNA候选序列
        步长: 1nt
        """
        if start_pos < 0:
            start_pos = 0
        if end_pos > len(target_mrna):
            end_pos = len(target_mrna)
        if start_pos >= end_pos or end_pos - start_pos < length:
            return []

        candidates = []
        mrna_length = len(target_mrna)
        
        target_mrna = target_mrna.replace('T', 'U').replace('t', 'u')

        for i in range(start_pos, end_pos - length + 1):
            sense = target_mrna[i:i+length].upper()

            if not self.is_valid_sequence(sense):
                continue

            evaluation = self.comprehensive_sirna_evaluation(sense, i, mrna_length)

            if not evaluation['passes_hard_filters']:
                continue

            evaluation['position'] = i + 1
            candidates.append(evaluation)

        candidates.sort(key=lambda x: x['final_score'], reverse=True)

        return candidates[:max_candidates]

    def check_sequence_quality(self, sequence: str) -> Dict[str, float]:
        results = {}
        gc = self.calculate_gc_content(sequence)
        results['gc_content'] = gc
        has_run = bool(re.search(r'GGGG|UUUU', sequence))
        results['has_run'] = has_run
        tm = SequenceAnalyzer.calculate_tm(sequence)
        results['tm'] = tm
        return results

    def get_specificity_score(self, sense: str, mrna: str) -> float:
        off_target_result = self.check_off_target(sense, mrna)
        if not off_target_result['has_off_target']:
            return 100.0
        penalty = min(len(off_target_result['matches']) * 5, 50)
        return max(0.0, 100.0 - penalty)

    def design_aso(self, target_mrna: str, target_region: str) -> Optional[str]:
        if target_region not in target_mrna:
            return None
        return self.reverse_complement(target_region)

    def optimize_sequence(self, sequence: str, target_gc: float = 50.0) -> str:
        current_gc = self.calculate_gc_content(sequence)
        optimized = list(sequence.upper())

        if current_gc < target_gc:
            for i in range(len(optimized)):
                if optimized[i] == 'A':
                    optimized[i] = 'G'
                elif optimized[i] == 'T':
                    optimized[i] = 'C'
                if self.calculate_gc_content(''.join(optimized)) >= target_gc:
                    break
        elif current_gc > target_gc:
            for i in range(len(optimized)):
                if optimized[i] == 'G':
                    optimized[i] = 'A'
                elif optimized[i] == 'C':
                    optimized[i] = 'T'
                if self.calculate_gc_content(''.join(optimized)) <= target_gc:
                    break

        return ''.join(optimized)

    def generate_shrna(self, sense_strand: str, loop_sequence: str = 'UUCG') -> str:
        antisense_strand = self.reverse_complement(sense_strand)
        return sense_strand + loop_sequence + antisense_strand


class SequenceAnalyzer:
    """
    序列分析工具
    """

    @staticmethod
    def find_motifs(sequence: str, motif: str) -> List[int]:
        positions = []
        motif_len = len(motif)
        for i in range(len(sequence) - motif_len + 1):
            if sequence[i:i+motif_len].upper() == motif.upper():
                positions.append(i+1)
        return positions

    @staticmethod
    def calculate_tm(sequence: str) -> float:
        tm = 64.9 + 41 * (sequence.count('G') + sequence.count('C') - 16.4) / len(sequence)
        return tm

    @staticmethod
    def detect_palindrome(sequence: str, min_length: int = 4) -> List[str]:
        palindromes = []
        for length in range(min_length, len(sequence)//2 + 1):
            for i in range(len(sequence) - length * 2 + 1):
                segment = sequence[i:i+length]
                complement_segment = NucleicAcidDesigner().reverse_complement(segment)
                if sequence[i+length:i+length*2] == complement_segment:
                    palindromes.append(sequence[i:i+length*2])
        return palindromes
