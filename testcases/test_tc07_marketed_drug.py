# -*- coding: utf-8 -*-
"""TC-07 全功能测试用例：已上市 siRNA 药物 Patisiran（Onpattro）

药物序列来自 IUPHAR/BPS Guide to Pharmacology 与 FDA 申报公开信息：
- Patisiran（ALN-18328），FDA 2018 批准，首个上市 siRNA 药物
- 靶点：TTR mRNA 3'UTR（突变型与野生型）；递送：LNP 肝靶向；适应症：hATTR
- 双链各 19nt RNA + 3' dTdT 悬垂，正义链 pos7/17 为 2'-OMe
用例文档见 TC-07-全功能测试-已上市药物Patisiran.md
"""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.chemical_modification import ChemicalModifier, DeliverySystem
from src.data_processing import DataHandler
from src.sequence_design import NucleicAcidDesigner, SequenceAnalyzer

# 药物档案（公开信息登记）
DRUG_PROFILE = {
    'generic_name': 'patisiran',
    'brand_name': 'Onpattro',
    'code': 'ALN-18328',
    'target': 'TTR mRNA 3UTR',
    'delivery': 'LNP',
    'indication': 'hATTR',
    'fda_approved_year': 2018,
}

SENSE = 'AUGGAAUACUCUUGGUUAC'      # passenger 链 5'->3'（19nt）
ANTISENSE = 'GUAACCAAGAGUAUUCCAU'  # guide 链 5'->3'（19nt）


class TestT07PatisiranFullFunction(unittest.TestCase):
    """以已上市药物真实序列串起全功能流程"""

    def setUp(self):
        self.d = NucleicAcidDesigner()

    def test_tc07_001_drug_profile_and_sequence(self):
        """TC07-001 药物档案与序列完整性"""
        self.assertEqual(DRUG_PROFILE['brand_name'], 'Onpattro')
        self.assertEqual(DRUG_PROFILE['delivery'], 'LNP')
        for seq in (SENSE, ANTISENSE):
            self.assertEqual(len(seq), 19)
            self.assertTrue(self.d.is_valid_sequence(seq))
            self.assertTrue(set(seq) <= {'A', 'U', 'G', 'C'})

    def test_tc07_002_duplex_complementarity(self):
        """TC07-002 双链互补性验证：guide = passenger 反向互补"""
        self.assertEqual(self.d.reverse_complement(SENSE), ANTISENSE)

    def test_tc07_003_physicochemical_properties(self):
        """TC07-003 理化性质评估"""
        gc = self.d.calculate_gc_content(SENSE)
        self.assertAlmostEqual(gc, 36.84, places=2)
        # GC 处于硬过滤 HF1 的 25-65% 区间
        self.assertTrue(25.0 <= gc <= 65.0)
        tm = SequenceAnalyzer.calculate_tm(SENSE)
        self.assertAlmostEqual(tm, 44.62, places=2)
        safe, warnings = self.d.check_immune_stimulatory(SENSE)
        self.assertTrue(safe)
        self.assertEqual(warnings, [])

    def test_tc07_004_comprehensive_evaluation(self):
        """TC07-004 综合评估（行为登记：真实裸序列得'良好'，热力学不对称判据不通过）"""
        ev = self.d.comprehensive_sirna_evaluation(SENSE)
        self.assertAlmostEqual(ev['final_score'], 69.52, places=2)
        self.assertEqual(ev['recommendation'], '良好')
        self.assertEqual(ev['reynolds_score'], 4)
        self.assertEqual(ev['amarzguioui_score'], 5)
        self.assertTrue(ev['uitei_passed'])
        self.assertFalse(ev['thermo_asymmetry'])
        self.assertTrue(ev['passes_hard_filters'])
        self.assertEqual(ev['antisense'], ANTISENSE)

    def test_tc07_005_lnp_modification_scheme(self):
        """TC07-005 LNP 修饰方案（与真实药物 dTdT 悬垂一致）"""
        scheme = ChemicalModifier().design_siRNA_modification(SENSE, ANTISENSE, platform='LNP')
        self.assertEqual(scheme['platform'], 'LNP')
        self.assertEqual(scheme['duplex_config'], '19+19')
        self.assertTrue(scheme['validation']['valid'])
        self.assertTrue(scheme['sense_modified'].endswith('-dTdT'))
        self.assertTrue(scheme['antisense_modified'].startswith('p'))
        # LNP 平台正义链仅奇数位 2'-OMe
        self.assertEqual(scheme['sense_map'][0]['sugar'], '2OMe')
        self.assertIsNone(scheme['sense_map'][1]['sugar'])

    def test_tc07_006_delivery_and_dosage(self):
        """TC07-006 递送方式与剂量（Patisiran 即 LNP 肝靶向）"""
        methods = DeliverySystem.recommend_delivery('肝脏')
        self.assertIn('LNP', methods)
        self.assertAlmostEqual(DeliverySystem.calculate_dosage(19, '肝脏'), 0.95, places=6)

    def test_tc07_007_end_to_end_archive(self):
        """TC07-007 全流程端到端存档"""
        # 评估 + 修饰方案
        ev = self.d.comprehensive_sirna_evaluation(SENSE)
        scheme = ChemicalModifier().design_siRNA_modification(SENSE, ANTISENSE, platform='LNP')
        # 药物档案存 JSON
        with tempfile.TemporaryDirectory() as tmp:
            profile_path = os.path.join(tmp, 'patisiran_profile.json')
            DataHandler.write_json(DRUG_PROFILE, profile_path)
            profile_loaded = DataHandler.read_json(profile_path)
            # 设计结果存 CSV
            results_path = os.path.join(tmp, 'patisiran_results.csv')
            DataHandler.save_design_results([{
                'drug': DRUG_PROFILE['generic_name'],
                'sense': SENSE,
                'antisense': ANTISENSE,
                'final_score': str(ev['final_score']),
                'recommendation': ev['recommendation'],
                'sense_modified': scheme['sense_modified'],
            }], results_path)
            rows = DataHandler.read_csv(results_path)
        # 读回校验
        self.assertEqual(profile_loaded['code'], 'ALN-18328')
        self.assertEqual(profile_loaded['fda_approved_year'], 2018)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['drug'], 'patisiran')
        self.assertEqual(rows[0]['sense'], SENSE)
        self.assertEqual(rows[0]['recommendation'], '良好')
        self.assertTrue(rows[0]['sense_modified'].endswith('-dTdT'))


if __name__ == '__main__':
    unittest.main(verbosity=2)
