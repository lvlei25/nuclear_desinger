# -*- coding: utf-8 -*-
"""TC-09 配置外置测试用例（ScoringConfig 默认值回归 + 自定义配置生效，用例文档见 TC-09-配置外置.md）"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config.settings import ScoringConfig
from src.sequence_design import NucleicAcidDesigner

# 与历史硬编码一致的默认配置基准
S1 = 'AUGAGGAGAUCUCCUCAGGCA'  # 21nt


class TestT09ScoringConfigDefaults(unittest.TestCase):
    """默认配置与历史硬编码值一致，保证行为不变"""

    def test_tc09_001_default_config_values(self):
        """TC09-001 默认配置字段值回归"""
        cfg = ScoringConfig.default()
        # 硬过滤
        self.assertEqual(cfg.hf1_gc_min, 25.0)
        self.assertEqual(cfg.hf1_gc_max, 65.0)
        self.assertEqual(cfg.hf2_poly_run, 5)
        self.assertEqual(cfg.hf3_reynolds_min, 2)
        # 脱靶
        self.assertEqual(cfg.off_target_identity_min, 80.0)
        self.assertEqual(cfg.off_target_identity_max, 100.0)
        self.assertEqual(cfg.off_target_penalty_per_match, 5.0)
        self.assertEqual(cfg.off_target_penalty_cap, 50.0)
        # 评分公共阈值与权重
        self.assertEqual(cfg.rule_gc_min, 30.0)
        self.assertEqual(cfg.rule_gc_max, 55.0)
        self.assertEqual(cfg.uitei_dg_threshold, -5.0)
        self.assertEqual(cfg.hard_pass_bonus, 30.0)
        self.assertEqual(cfg.reynolds_weight, 4.0)
        self.assertEqual(cfg.amarzguioui_weight, 3.0)
        self.assertEqual(cfg.seed_complexity_weight, 8.0)
        self.assertEqual(cfg.thermo_asymmetry_bonus, 3.0)
        self.assertEqual(cfg.uitei_pass_bonus, 3.0)
        self.assertEqual(cfg.optimal_gc_bonus, 5.0)
        self.assertEqual(cfg.max_score, 96.0)
        self.assertEqual(cfg.recommend_threshold, 70.0)
        self.assertEqual(cfg.good_threshold, 50.0)
        # 靶区域
        self.assertEqual(cfg.utr3_bonus, 10.0)
        self.assertEqual(cfg.cds_bonus, 5.0)
        self.assertEqual(cfg.utr5_bonus, 0.0)
        self.assertEqual(cfg.start_codon_penalty, -5.0)
        self.assertEqual(cfg.cds_start_fraction, 0.25)
        self.assertEqual(cfg.cds_end_fraction, 0.75)
        self.assertEqual(cfg.start_codon_window, 50)
        # 比对
        self.assertEqual(cfg.seed_size, 7)
        self.assertEqual(cfg.seed_flank, 5)

    def test_tc09_002_designer_uses_default_config(self):
        """TC09-002 设计器默认持有默认配置实例"""
        d = NucleicAcidDesigner()
        self.assertIsInstance(d.config, ScoringConfig)
        self.assertEqual(d.config.hf1_gc_min, 25.0)
        # 默认构造等价于显式传默认配置
        d2 = NucleicAcidDesigner(config=ScoringConfig.default())
        self.assertEqual(
            d.comprehensive_sirna_evaluation(S1)['final_score'],
            d2.comprehensive_sirna_evaluation(S1)['final_score'],
        )


class TestT09CustomConfig(unittest.TestCase):
    """自定义配置生效，默认配置不受影响"""

    def test_tc09_003_hard_filter_uses_custom_gc_range(self):
        """TC09-003 自定义 HF1 GC 区间生效"""
        d_default = NucleicAcidDesigner()
        # S1 GC = 10/21 ≈ 47.6%，落在默认 25-65 与收紧后的 48-65 判定差异之外/之内的构造不便，
        # 改用极宽区间 [0, 100] 与极窄区间验证判定翻转
        gc = d_default.calculate_gc_content(S1)

        wide = NucleicAcidDesigner(config=ScoringConfig(hf1_gc_min=0.0, hf1_gc_max=100.0))
        narrow = NucleicAcidDesigner(config=ScoringConfig(hf1_gc_min=gc + 1, hf1_gc_max=gc + 2))

        all_pass_wide, res_wide = wide.hard_filter(S1)
        _, res_narrow = narrow.hard_filter(S1)
        # 宽区间下 HF1 必过；窄区间下 HF1 必不过
        self.assertTrue(res_wide['hf1_gc']['pass'])
        self.assertFalse(res_narrow['hf1_gc']['pass'])
        # 标准文案跟随配置
        self.assertEqual(res_wide['hf1_gc']['standard'], '0-100%')

    def test_tc09_004_off_target_threshold_configurable(self):
        """TC09-004 自定义脱靶判定阈值生效"""
        d = NucleicAcidDesigner()
        # 结果中的 threshold 字段反映当前配置
        result = d.check_off_target(S1, S1)
        self.assertEqual(result['threshold'], 80.0)

        custom = NucleicAcidDesigner(config=ScoringConfig(off_target_identity_min=60.0))
        result_custom = custom.check_off_target(S1, S1)
        self.assertEqual(result_custom['threshold'], 60.0)

    def test_tc09_005_recommendation_thresholds_configurable(self):
        """TC09-005 自定义推荐分级阈值生效"""
        d_default = NucleicAcidDesigner()
        score = d_default.comprehensive_sirna_evaluation(S1)['final_score']

        # 把"推荐"下限抬到满分以上，应降级为"良好"或"一般"
        high_bar = NucleicAcidDesigner(config=ScoringConfig(recommend_threshold=score + 1))
        rec = high_bar.comprehensive_sirna_evaluation(S1)['recommendation']
        self.assertIn(rec, ('良好', '一般'))

        # 把分级下限归零，硬过滤通过时应为"推荐"
        low_bar = NucleicAcidDesigner(config=ScoringConfig(recommend_threshold=0.0, good_threshold=0.0))
        self.assertEqual(low_bar.comprehensive_sirna_evaluation(S1)['recommendation'], '推荐')

    def test_tc09_006_config_isolation_between_instances(self):
        """TC09-006 不同实例配置相互隔离"""
        d1 = NucleicAcidDesigner(config=ScoringConfig(seed_size=5))
        d2 = NucleicAcidDesigner()
        self.assertEqual(d1.config.seed_size, 5)
        self.assertEqual(d2.config.seed_size, 7)


if __name__ == '__main__':
    unittest.main(verbosity=2)
