#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, 'c:/Users/lei.lv/ide/小核酸药物设计')

from src.sequence_design import NucleicAcidDesigner

# 创建设计器实例
designer = NucleicAcidDesigner()

# 用户提供的测试序列
sense_sequence = "UCACGACCGUGUACAUUGACC"
antisense_sequence = "GGUCAAUGUACACGGUCGUGA"
target_region = "3'-UTR"

print("=" * 70)
print("序列评分验证测试")
print("=" * 70)
print("正义链:", sense_sequence)
print("反义链:", antisense_sequence)
print("靶区域:", target_region)
print("-" * 70)

# 调用评价方法
evaluation = designer.comprehensive_sirna_evaluation(
    sense_sequence, 
    target_region_name=target_region
)

print("评价结果:")
print("GC含量: %.1f%%" % evaluation['gc_content'])
print("Reynolds评分: %d/8" % evaluation['reynolds_score'])
print("Amarzguioui评分: %d/5" % evaluation['amarzguioui_score'])
print("种子区复杂度: %.2f" % evaluation['seed_complexity'])
print("热力学不对称性: %s" % ("通过" if evaluation['thermo_asymmetry'] else "不通过"))
print("Ui-Tei规则: %s" % ("通过" if evaluation['uitei_passed'] else "不通过"))
print("靶区域: %s" % evaluation['target_region'])
print("靶区域加分: %d" % evaluation['target_bonus'])
print("最终总分: %.1f/%d" % (evaluation['final_score'], evaluation['max_score']))
print("推荐等级: %s" % evaluation['recommendation'])
print("-" * 70)

# 详细计算过程
print("评分计算明细:")
print("  1. 硬性过滤通过: 30分")
print("  2. Reynolds评分(%d × 4): %d分" % (evaluation['reynolds_score'], evaluation['reynolds_score'] * 4))
print("  3. Amarzguioui评分(%d × 3): %d分" % (evaluation['amarzguioui_score'], evaluation['amarzguioui_score'] * 3))
print("  4. 种子区复杂度(%.2f × 8): %.2f分" % (evaluation['seed_complexity'], evaluation['seed_complexity'] * 8))
print("  5. 热力学不对称性(通过+3): %d分" % (3 if evaluation['thermo_asymmetry'] else 0))
print("  6. Ui-Tei规则(通过+3): %d分" % (3 if evaluation['uitei_passed'] else 0))
print("  7. GC含量最优(通过+5): %d分" % (5 if evaluation['optimal_gc'] else 0))
print("  8. 靶区域加分(%s): %d分" % (evaluation['target_region'], evaluation['target_bonus']))
print("-" * 70)
total = (30 if evaluation['passes_hard_filters'] else 0) + \
        evaluation['reynolds_score'] * 4 + \
        evaluation['amarzguioui_score'] * 3 + \
        evaluation['seed_complexity'] * 8 + \
        (3 if evaluation['thermo_asymmetry'] else 0) + \
        (3 if evaluation['uitei_passed'] else 0) + \
        (5 if evaluation['optimal_gc'] else 0) + \
        evaluation['target_bonus']
print("总分计算: %.2f" % total)
print("=" * 70)

# 验证是否与用户提供的结果一致
expected_score = 90.8
actual_score = evaluation['final_score']
if abs(actual_score - expected_score) < 0.1:
    print("评分验证通过！计算结果(%.1f)与预期(%.1f)一致" % (actual_score, expected_score))
else:
    print("评分验证失败！计算结果(%.1f)与预期(%.1f)不一致" % (actual_score, expected_score))