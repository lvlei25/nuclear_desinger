# -*- coding: utf-8 -*-
"""TC-05 综合评估与 GUI 集成测试用例（用例文档见 TC-05-综合评估与GUI.md）"""

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.chemical_modification import ChemicalModifier
from src.data_processing import DataHandler
from src.sequence_design import NucleicAcidDesigner

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S1 = 'AUAGCUGCAUGAAAGUGCGCC'
S3 = 'UUUUUUGGGGGCCCCCCCCCA'


class TestT05ComprehensiveEvaluation(unittest.TestCase):
    """综合评估"""

    def setUp(self):
        self.d = NucleicAcidDesigner()

    def test_tc05_001_evaluation_default_region(self):
        """TC05-001 综合评估-默认区域"""
        ev = self.d.comprehensive_sirna_evaluation(S1)
        self.assertAlmostEqual(ev['final_score'], 76.37, places=2)
        self.assertEqual(ev['recommendation'], '推荐')
        self.assertEqual(ev['target_region'], '起始密码子附近')
        self.assertEqual(ev['target_bonus'], -5)

    def test_tc05_002_evaluation_3utr(self):
        """TC05-002 综合评估-指定3'-UTR"""
        ev = self.d.comprehensive_sirna_evaluation(S1, target_region_name="3'-UTR")
        self.assertAlmostEqual(ev['final_score'], 91.37, places=2)
        self.assertEqual(ev['target_bonus'], 10)
        self.assertEqual(ev['recommendation'], '推荐')

    def test_tc05_003_evaluation_not_recommended(self):
        """TC05-003 综合评估-硬过滤失败不推荐"""
        ev = self.d.comprehensive_sirna_evaluation(S3)
        self.assertFalse(ev['passes_hard_filters'])
        self.assertEqual(ev['recommendation'], '不推荐')

    def test_tc05_004_antisense_consistency(self):
        """TC05-004 反义链与反向互补一致性"""
        ev = self.d.comprehensive_sirna_evaluation(S1)
        self.assertEqual(ev['antisense'], self.d.reverse_complement(S1))
        self.assertEqual(ev['antisense'], 'GGCGCACUUUCAUGCAGCUAU')


class TestT05Smoke(unittest.TestCase):
    """demo 与 GUI 冒烟"""

    def test_tc05_005_demo_smoke(self):
        """TC05-005 demo命令行冒烟"""
        result = subprocess.run(
            [sys.executable, 'demo.py'],
            cwd=ROOT, capture_output=True, text=True, timeout=120,
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr[-500:] if result.stderr else '')

    def test_tc05_006_gui_import_smoke(self):
        """TC05-006 GUI模块导入冒烟（不启动主窗口）"""
        result = subprocess.run(
            [sys.executable, '-c', 'import gui; print("GUI_IMPORT_OK")'],
            cwd=ROOT, capture_output=True, text=True, timeout=60,
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr[-500:] if result.stderr else '')
        self.assertIn('GUI_IMPORT_OK', result.stdout)


class TestT05EndToEnd(unittest.TestCase):
    """端到端流程：设计 -> 评估 -> 修饰 -> 存取"""

    def test_tc05_007_end_to_end(self):
        """TC05-007 端到端流程"""
        # 1. 设计 siRNA
        d = NucleicAcidDesigner()
        mrna = S1 * 5
        sense, antisense = d.design_sirna(mrna, position=30, length=21)
        # 2. 综合评估
        ev = d.comprehensive_sirna_evaluation(sense)
        self.assertIn(ev['recommendation'], ('推荐', '良好', '一般'))
        # 3. GalNAc 修饰方案
        scheme = ChemicalModifier().design_siRNA_modification(sense, antisense)
        self.assertTrue(scheme['validation']['valid'])
        # 4. 结果存 CSV 并读回
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'e2e_results.csv')
            DataHandler.save_design_results(
                [{'sense': sense, 'final_score': str(ev['final_score']),
                  'sense_modified': scheme['sense_modified']}],
                path,
            )
            with open(path, 'r', encoding='utf-8') as f:
                rows = list(csv.DictReader(f))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['sense'], sense)
        self.assertIn('final_score', rows[0])
        self.assertIn('sense_modified', rows[0])


def _url_ctx(payload):
    """构造 Mock 的 urlopen 上下文管理器返回值（JSON 字典或字节文本）"""
    resp = MagicMock()
    resp.__enter__.return_value = resp
    if isinstance(payload, bytes):
        resp.read.return_value = payload
    else:
        resp.read.return_value = json.dumps(payload).encode('utf-8')
    return resp


class TestT05NCBIAlias(unittest.TestCase):
    """基因别名自动解析（gui.PubMedQuerier，requests/urllib 全部 Mock，不联网）"""

    def test_tc05_008_resolve_official_symbol(self):
        """TC05-008 别名解析单元：Gene库 ESearch -> ESummary 取官方符号"""
        from gui import PubMedQuerier
        payloads = [
            {'esearchresult': {'idlist': ['630504']}},
            {'result': {'630504': {'name': 'ACVR1C'}}},
        ]
        with patch('gui.urllib.request.urlopen', side_effect=[_url_ctx(p) for p in payloads]):
            symbol = PubMedQuerier.resolve_official_symbol('ALK-7', 'Homo sapiens[Organism]')
        self.assertEqual(symbol, 'ACVR1C')

    def test_tc05_009_search_fallback_with_alias(self):
        """TC05-009 检索回退流程：直接检索无命中 -> 别名解析 -> 官方符号重试 -> 取序列"""
        from gui import PubMedQuerier
        fasta = b'>NM_1|Homo sapiens ACVR1C (activin A receptor type 1C), mRNA\nACGUACGU\n'
        payloads = [
            {'esearchresult': {'idlist': []}},                      # 1. Nucleotide 直接检索无命中
            {'esearchresult': {'idlist': ['630504']}},              # 2. Gene 库解析别名
            {'result': {'630504': {'name': 'ACVR1C'}}},             # 3. ESummary 取官方符号
            {'esearchresult': {'idlist': ['111', '222']}},          # 4. 官方符号重试检索命中
            fasta,                                                  # 5. EFetch 返回 FASTA 文本
        ]
        with patch('gui.urllib.request.urlopen', side_effect=[_url_ctx(p) for p in payloads]) as mock_urlopen:
            sequences = PubMedQuerier.search_nucleotide_mrna('ALK-7', max_results=5)
        self.assertEqual(len(sequences), 1)
        self.assertEqual(sequences[0]['accession'], 'NM_1')
        self.assertIn('ACVR1C', sequences[0]['definition'])
        self.assertEqual(sequences[0]['sequence'], 'ACGUACGU')
        # 第4次请求应使用官方符号 ACVR1C 检索
        self.assertIn('ACVR1C', mock_urlopen.call_args_list[3][0][0])


class TestT05GuiSourceAudit(unittest.TestCase):
    """GUI 源码静态审计（防回归）"""

    def test_tc05_010_no_undefined_gene_var(self):
        """TC05-010 gui.py 不存在未定义的 gene_name_var 引用"""
        # 2026-10-15 缺陷：2660/2755 引用从未定义的 self.gene_name_var，导致 siRNA 规则评价中断
        src_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'gui.py')
        with open(src_path, 'r', encoding='utf-8') as f:
            source = f.read()
        self.assertNotIn('self.gene_name_var', source)
        self.assertIn('self.gene_id_entry.get().strip()', source)


if __name__ == '__main__':
    unittest.main(verbosity=2)
