# -*- coding: utf-8 -*-
"""TC-03 数据处理模块测试用例（用例文档见 TC-03-数据处理.md）"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.data_processing import DataHandler, SequenceDatabase


class TestT03DataHandler(unittest.TestCase):
    """FASTA/CSV/JSON 读写"""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.addCleanup(self._tmp.cleanup)

    def test_tc03_001_fasta_roundtrip(self):
        """TC03-001 FASTA写读回环与60nt换行"""
        seq_long = 'ACGT' * 33  # 132nt
        data = {'seq_short': 'AUGCAU', 'seq_long': seq_long}
        path = os.path.join(self.tmp, 't.fasta')
        DataHandler.write_fasta(data, path)
        with open(path, 'r', encoding='utf-8') as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        # 132nt 拆为 60+60+12 三行；两条序列共 2 个头 + 1 + 3 = 6 行
        self.assertEqual(len(lines), 6)
        self.assertEqual(len(lines[3]), 60)
        self.assertEqual(len(lines[5]), 12)
        loaded = DataHandler.read_fasta(path)
        self.assertEqual(loaded, data)

    def test_tc03_002_fasta_empty_header_ignored(self):
        """TC03-002 FASTA空序列头忽略"""
        path = os.path.join(self.tmp, 'empty.fasta')
        with open(path, 'w', encoding='utf-8') as f:
            f.write('>empty_head\n>real_head\nAUGCAU\n')
        loaded = DataHandler.read_fasta(path)
        self.assertEqual(loaded, {'real_head': 'AUGCAU'})

    def test_tc03_003_csv_roundtrip(self):
        """TC03-003 CSV写读回环（含中文）"""
        data = [
            {'name': '候选一', 'sense': 'AUGCAU', 'score': '88.5'},
            {'name': '候选二', 'sense': 'GGCCAA', 'score': '76.0'},
        ]
        path = os.path.join(self.tmp, 't.csv')
        DataHandler.write_csv(data, path)
        loaded = DataHandler.read_csv(path)
        self.assertEqual(loaded, data)

    def test_tc03_004_csv_empty_no_file(self):
        """TC03-004 CSV空数据不生成文件"""
        path = os.path.join(self.tmp, 'empty.csv')
        DataHandler.write_csv([], path)
        self.assertFalse(os.path.exists(path))

    def test_tc03_005_json_roundtrip(self):
        """TC03-005 JSON写读回环（嵌套与中文）"""
        data = {'name': '测试序列', 'tags': ['肝靶向', '稳定'], 'meta': {'length': 21, 'gc': 52.38}}
        path = os.path.join(self.tmp, 't.json')
        DataHandler.write_json(data, path)
        with open(path, 'r', encoding='utf-8') as f:
            raw = f.read()
        self.assertIn('测试序列', raw)  # ensure_ascii=False，中文不转义
        self.assertEqual(DataHandler.read_json(path), data)

    def test_tc03_006_load_target_genes(self):
        """TC03-006 目标基因加载别名"""
        data = [{'gene': 'TTR', 'mrna': 'AUGCAU'}]
        path = os.path.join(self.tmp, 'targets.csv')
        DataHandler.write_csv(data, path)
        self.assertEqual(DataHandler.load_target_genes(path), data)

    def test_tc03_007_save_design_results(self):
        """TC03-007 设计结果保存别名"""
        data = [{'sense': 'AUGCAU', 'final_score': '76.37'}]
        path = os.path.join(self.tmp, 'results.csv')
        DataHandler.save_design_results(data, path)
        self.assertEqual(DataHandler.read_csv(path), data)


class TestT03SequenceDatabase(unittest.TestCase):
    """序列数据库"""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.db = SequenceDatabase(db_path=os.path.join(self._tmp.name, 'db'))
        self.addCleanup(self._tmp.cleanup)

    def test_tc03_008_database_crud(self):
        """TC03-008 序列库增查列删回环"""
        # 增
        self.db.add_sequence('test_seq', 'AUGCAUGC', '测试描述')
        # 查
        got = self.db.get_sequence('test_seq')
        self.assertIsNotNone(got)
        self.assertEqual(got['name'], 'test_seq')
        self.assertEqual(got['sequence'], 'AUGCAUGC')
        self.assertEqual(got['length'], 8)
        self.assertEqual(got['description'], '测试描述')
        # 列
        self.assertIn('test_seq', self.db.list_sequences())
        # 删
        self.assertTrue(self.db.delete_sequence('test_seq'))
        self.assertIsNone(self.db.get_sequence('test_seq'))
        self.assertFalse(self.db.delete_sequence('test_seq'))
        # 查询不存在
        self.assertIsNone(self.db.get_sequence('not_exist'))


if __name__ == '__main__':
    unittest.main(verbosity=2)
