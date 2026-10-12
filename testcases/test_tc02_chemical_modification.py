# -*- coding: utf-8 -*-
"""TC-02 化学修饰模块测试用例（用例文档见 TC-02-化学修饰.md）"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.chemical_modification import ChemicalModifier, DeliverySystem
from src.sequence_design import NucleicAcidDesigner

S1 = 'AUAGCUGCAUGAAAGUGCGCC'


class TestT02ModificationOps(unittest.TestCase):
    """修饰的增删改查与输出格式"""

    def setUp(self):
        self.cm = ChemicalModifier()

    def test_tc02_001_add_modification(self):
        """TC02-001 添加修饰"""
        self.cm.add_modification(2, '2OMe')
        mods = self.cm.get_modifications()
        self.assertEqual(len(mods), 1)
        self.assertEqual(mods[0]['position'], 2)
        self.assertEqual(mods[0]['type'], '2OMe')
        self.assertIn('effect', mods[0]['info'])

    def test_tc02_002_unknown_mod_type(self):
        """TC02-002 未知修饰类型抛ValueError"""
        with self.assertRaises(ValueError):
            self.cm.add_modification(1, 'XYZ')

    def test_tc02_003_remove_and_clear(self):
        """TC02-003 移除与清空修饰"""
        self.cm.add_modification(1, 'PS')
        self.cm.add_modification(2, '2OMe')
        self.cm.remove_modification(1)
        self.assertEqual(len(self.cm.get_modifications()), 1)
        self.cm.clear_modifications()
        self.assertEqual(self.cm.get_modifications(), [])

    def test_tc02_004_apply_mark_format(self):
        """TC02-004 标记格式输出"""
        self.cm.add_modification(2, '2OMe')
        self.assertEqual(self.cm.apply_modifications('ATCG'), 'ATmCG')

    def test_tc02_005_apply_synth_format(self):
        """TC02-005 合成订单格式输出"""
        self.cm.add_modification(2, '2OMe')
        self.assertEqual(self.cm.apply_modifications('ATCG', 'synth'), 'AmTCG')

    def test_tc02_006_apply_out_of_range(self):
        """TC02-006 越界位置不生效"""
        self.cm.add_modification(10, '2OMe')
        self.assertEqual(self.cm.apply_modifications('ATCG'), 'ATCG')


class TestT02Suggestion(unittest.TestCase):
    """修饰策略推荐"""

    def test_tc02_007_suggest_sirna(self):
        """TC02-007 siRNA推荐策略"""
        suggestions = ChemicalModifier().suggest_modifications('siRNA')
        self.assertEqual(len(suggestions), 4)
        types = {s['type'] for s in suggestions}
        self.assertIn('2OMe', types)
        self.assertIn('PS', types)

    def test_tc02_008_suggest_aso_mirna(self):
        """TC02-008 ASO/miRNA推荐策略"""
        self.assertGreater(len(ChemicalModifier().suggest_modifications('ASO')), 0)
        self.assertGreater(len(ChemicalModifier().suggest_modifications('miRNA')), 0)


class TestT02DesignScheme(unittest.TestCase):
    """siRNA 修饰方案设计（GalNAc / LNP）"""

    @classmethod
    def setUpClass(cls):
        cls.sense = S1
        cls.antisense = NucleicAcidDesigner().reverse_complement(S1)

    def test_tc02_009_galnac_scheme(self):
        """TC02-009 GalNAc修饰方案"""
        result = ChemicalModifier().design_siRNA_modification(self.sense, self.antisense)
        self.assertEqual(result['platform'], 'GalNAc')
        self.assertEqual(result['duplex_config'], '21+21')
        self.assertTrue(result['sense_modified'].endswith('-GalNAc'))
        self.assertTrue(result['antisense_modified'].startswith('p'))
        self.assertTrue(result['validation']['valid'])
        # GalNAc 平台正义链全链 2'-OMe
        self.assertTrue(all(m['sugar'] == '2OMe' for m in result['sense_map']))

    def test_tc02_010_lnp_scheme(self):
        """TC02-010 LNP修饰方案"""
        result = ChemicalModifier().design_siRNA_modification(self.sense, self.antisense, platform='LNP')
        self.assertEqual(result['platform'], 'LNP')
        self.assertTrue(result['sense_modified'].endswith('-dTdT'))
        # LNP 平台仅奇数位覆盖 2'-OMe
        self.assertEqual(result['sense_map'][0]['sugar'], '2OMe')
        self.assertIsNone(result['sense_map'][1]['sugar'])
        self.assertEqual(result['sense_map'][2]['sugar'], '2OMe')

    def test_tc02_011_seed_coverage(self):
        """TC02-011 种子区糖修饰覆盖（AS pos2-8）"""
        result = ChemicalModifier().design_siRNA_modification(self.sense, self.antisense)
        seed_sugars = [result['antisense_map'][i]['sugar'] for i in range(1, 8)]
        self.assertTrue(all(s is not None for s in seed_sugars))

    def test_tc02_012_cut_center_no_fana(self):
        """TC02-012 切割位点中心区排除 FANA（AS 索引 8-11）"""
        result = ChemicalModifier().design_siRNA_modification(self.sense, self.antisense)
        center = [result['antisense_map'][i]['sugar'] for i in range(8, 12)]
        self.assertNotIn('FANA', center)


class TestT02Delivery(unittest.TestCase):
    """递送系统与剂量"""

    def test_tc02_013_recommend_delivery(self):
        """TC02-013 递送方式推荐"""
        self.assertEqual(DeliverySystem.recommend_delivery('肝脏'), ['GalNAc', 'LNP'])
        self.assertEqual(DeliverySystem.recommend_delivery('未知组织'), ['LNP'])

    def test_tc02_014_calculate_dosage(self):
        """TC02-014 剂量计算"""
        self.assertAlmostEqual(DeliverySystem.calculate_dosage(21, '肝脏'), 1.05, places=6)
        self.assertAlmostEqual(DeliverySystem.calculate_dosage(40, '眼睛'), 1.0, places=6)


if __name__ == '__main__':
    unittest.main(verbosity=2)
