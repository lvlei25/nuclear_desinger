# -*- coding: utf-8 -*-
"""TC-01 序列设计模块测试用例（用例文档见 TC-01-序列设计.md）"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.sequence_design import NucleicAcidDesigner, SequenceAnalyzer

# 公共测试序列
S1 = 'AUAGCUGCAUGAAAGUGCGCC'   # 21nt，GC 52.38%，通过全部硬过滤
S2 = 'AUAGCUGCAUGCCAGUGCGAU'   # 21nt，Ui-Tei 反义链 5' 端不满足
S3 = 'UUUUUUGGGGGCCCCCCCCCA'   # 21nt，四条硬过滤全部不通过
MUT = 'AUAGCAUCAGCAAAGUGCGCC'  # S1 的突变版本（脱靶场景）


class TestT01BasicCalc(unittest.TestCase):
    """基础序列计算"""

    def setUp(self):
        self.d = NucleicAcidDesigner()

    def test_tc01_001_complement(self):
        """TC01-001 互补链计算"""
        self.assertEqual(self.d.complement('ATCG'), 'UAGC')

    def test_tc01_002_reverse_complement(self):
        """TC01-002 反向互补计算"""
        self.assertEqual(self.d.reverse_complement('ATCG'), 'CGAU')

    def test_tc01_003_reverse_complement_lowercase(self):
        """TC01-003 反向互补小写输入（行为登记：小写不转换）"""
        self.assertEqual(self.d.reverse_complement('atcg'), 'gcta')

    def test_tc01_004_gc_content(self):
        """TC01-004 GC含量计算"""
        self.assertAlmostEqual(self.d.calculate_gc_content('ATCGCG'), 66.67, places=2)

    def test_tc01_005_gc_content_empty(self):
        """TC01-005 GC含量空序列（P0-4 边界统一：抛 ValueError）"""
        with self.assertRaises(ValueError):
            self.d.calculate_gc_content('')

    def test_tc01_006_is_valid_sequence(self):
        """TC01-006 序列合法性校验（空串视为合法，行为登记）"""
        self.assertTrue(self.d.is_valid_sequence('AUCG'))
        self.assertFalse(self.d.is_valid_sequence('ATCGXYZ'))
        self.assertTrue(self.d.is_valid_sequence(''))


class TestT01Design(unittest.TestCase):
    """siRNA/ASO/shRNA/优化设计"""

    def setUp(self):
        self.d = NucleicAcidDesigner()

    def test_tc01_007_design_sirna_ok(self):
        """TC01-007 siRNA设计正常"""
        mrna = 'A' * 200
        sense, antisense = self.d.design_sirna(mrna, position=100, length=21)
        self.assertEqual(len(sense), 21)
        self.assertEqual(len(antisense), 21)
        self.assertEqual(antisense, self.d.reverse_complement(sense))

    def test_tc01_008_design_sirna_out_of_range(self):
        """TC01-008 siRNA设计越界抛ValueError"""
        with self.assertRaises(ValueError):
            self.d.design_sirna('A' * 200, position=190, length=21)

    def test_tc01_026_design_aso(self):
        """TC01-026 ASO设计"""
        self.assertEqual(self.d.design_aso(S1, S1[:8]), 'GCAGCUAU')
        self.assertIsNone(self.d.design_aso(S1, 'ZZZ'))

    def test_tc01_027_optimize_sequence(self):
        """TC01-027 序列GC优化"""
        self.assertEqual(self.d.optimize_sequence('ATAT'), 'GCAT')

    def test_tc01_028_generate_shrna(self):
        """TC01-028 shRNA生成"""
        self.assertEqual(self.d.generate_shrna('AUGCAAU'), 'AUGCAAUUUCGAUUGCAU')


class TestT01Scoring(unittest.TestCase):
    """四套评分规则与硬过滤"""

    def setUp(self):
        self.d = NucleicAcidDesigner()

    def test_tc01_009_reynolds_score(self):
        """TC01-009 Reynolds评分"""
        self.assertEqual(self.d.reynolds_score(S1), 6)
        self.assertEqual(self.d.reynolds_score(S2), 6)

    def test_tc01_010_amarzguioui_score(self):
        """TC01-010 Amarzguioui评分"""
        self.assertEqual(self.d.amarzguioui_score(S1), 3)

    def test_tc01_011_uitei_pass(self):
        """TC01-011 Ui-Tei规则通过"""
        passed, reason = self.d.uitei_check(S1)
        self.assertTrue(passed)
        self.assertEqual(reason, '完全符合Ui-Tei规则')

    def test_tc01_012_uitei_fail(self):
        """TC01-012 Ui-Tei规则失败原因"""
        passed, reason = self.d.uitei_check(S2)
        self.assertFalse(passed)
        self.assertIn('反义链1位为A', reason)

    def test_tc01_013_hard_filter_pass(self):
        """TC01-013 硬过滤全部通过"""
        all_pass, hf = self.d.hard_filter(S1)
        self.assertTrue(all_pass)
        self.assertTrue(all(hf[k]['pass'] for k in ('hf1_gc', 'hf2_poly', 'hf3_reynolds', 'hf4_immune')))

    def test_tc01_014_hard_filter_fail(self):
        """TC01-014 硬过滤全部失败"""
        all_pass, hf = self.d.hard_filter(S3)
        self.assertFalse(all_pass)
        self.assertFalse(hf['hf1_gc']['pass'])
        self.assertFalse(hf['hf2_poly']['pass'])
        self.assertFalse(hf['hf3_reynolds']['pass'])
        self.assertFalse(hf['hf4_immune']['pass'])

    def test_tc01_015_thermodynamic_asymmetry(self):
        """TC01-015 热力学不对称性"""
        passed, diff = self.d.thermodynamic_asymmetry(S1)
        self.assertTrue(passed)
        self.assertAlmostEqual(diff, -3.3, places=6)

    def test_tc01_016_seed_complexity(self):
        """TC01-016 种子区复杂度"""
        self.assertAlmostEqual(self.d.calculate_seed_complexity(S1), 0.9212, places=4)

    def test_tc01_017_immune_detected(self):
        """TC01-017 免疫刺激序列检出"""
        safe, warnings = self.d.check_immune_stimulatory('ACUGUGUAC')
        self.assertFalse(safe)
        self.assertTrue(any('UGUGU' in w for w in warnings))

    def test_tc01_018_immune_safe(self):
        """TC01-018 免疫刺激序列安全"""
        safe, warnings = self.d.check_immune_stimulatory(S1)
        self.assertTrue(safe)
        self.assertEqual(warnings, [])

    def test_tc01_031_short_sequence_scoring_guard(self):
        """TC01-031 短序列评分防护（P0-4：IndexError 统一为带明确信息的 ValueError）"""
        for seq in ('AUGC', ''):
            with self.assertRaises(ValueError) as ctx:
                self.d.reynolds_score(seq)
            self.assertIn('19nt', str(ctx.exception))
            with self.assertRaises(ValueError):
                self.d.amarzguioui_score(seq)


class TestT01RegionAndOffTarget(unittest.TestCase):
    """靶区域加成、脱靶与比对"""

    def setUp(self):
        self.d = NucleicAcidDesigner()

    def test_tc01_019_target_region_bonus_by_position(self):
        """TC01-019 靶区域加成-位置"""
        self.assertEqual(self.d.get_target_region_bonus(20, 400), (0, "5'-UTR"))
        self.assertEqual(self.d.get_target_region_bonus(100, 400), (-5, '起始密码子附近'))
        self.assertEqual(self.d.get_target_region_bonus(200, 400), (5, 'CDS'))
        self.assertEqual(self.d.get_target_region_bonus(350, 400), (10, "3'-UTR"))

    def test_tc01_020_target_region_bonus_by_name(self):
        """TC01-020 靶区域加成-名称"""
        self.assertEqual(self.d.get_target_region_bonus_by_name("3'-UTR"), (10, "3'-UTR"))
        self.assertEqual(self.d.get_target_region_bonus_by_name('未知区域'), (5, 'CDS'))

    def test_tc01_021_check_off_target(self):
        """TC01-021 脱靶检测"""
        # 完全匹配位点不计入脱靶（阈值区间为 [80, 100)）
        result_none = self.d.check_off_target(S1, S1 * 3)
        self.assertFalse(result_none['has_off_target'])
        # 突变链产生一处约 80.95% 相似度的脱靶
        result_hit = self.d.check_off_target(S1, S1 + MUT + S1)
        self.assertTrue(result_hit['has_off_target'])
        self.assertEqual(len(result_hit['matches']), 1)
        self.assertAlmostEqual(result_hit['max_identity'], 80.95, places=2)

    def test_tc01_022_specificity_score(self):
        """TC01-022 特异性评分"""
        self.assertAlmostEqual(self.d.get_specificity_score(S1, S1 + MUT + S1), 95.0, places=6)
        self.assertAlmostEqual(self.d.get_specificity_score(S1, S1 * 3), 100.0, places=6)

    def test_tc01_023_blast_alignment(self):
        """TC01-023 本地BLAST比对"""
        score, detail = self.d.blast_alignment(S1, S1)
        self.assertAlmostEqual(score, 100.0, places=6)
        score2, detail2 = self.d.blast_alignment('', 'AUGC')
        self.assertEqual(score2, 0.0)
        self.assertEqual(detail2, '序列为空')

    def test_tc01_024_find_sirna_candidates(self):
        """TC01-024 siRNA候选查找"""
        candidates = self.d.find_sirna_candidates(S1 * 10, 0, 210)
        self.assertEqual(len(candidates), 20)
        scores = [c['final_score'] for c in candidates]
        self.assertEqual(scores, sorted(scores, reverse=True))
        self.assertAlmostEqual(max(scores), 91.37, places=2)

    def test_tc01_025_find_sirna_candidates_invalid_range(self):
        """TC01-025 候选查找无效区间"""
        self.assertEqual(self.d.find_sirna_candidates(S1, 10, 5), [])
        self.assertEqual(self.d.find_sirna_candidates(S1, 0, 10), [])


class TestT01Analyzer(unittest.TestCase):
    """SequenceAnalyzer 与质量检查"""

    def setUp(self):
        self.d = NucleicAcidDesigner()

    def test_tc01_029_check_sequence_quality(self):
        """TC01-029 序列质量检查"""
        quality = self.d.check_sequence_quality(S1)
        self.assertAlmostEqual(quality['gc_content'], 52.38, places=2)
        self.assertFalse(quality['has_run'])
        self.assertAlmostEqual(quality['tm'], 54.36, places=2)

    def test_tc01_030_calculate_tm(self):
        """TC01-030 Tm值计算"""
        self.assertAlmostEqual(SequenceAnalyzer.calculate_tm('ATCGCGATCG'), 22.26, places=2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
