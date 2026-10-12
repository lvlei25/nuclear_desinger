# -*- coding: utf-8 -*-
"""TC-06 NCBI 在线访问集成测试用例（真实联网，网络不可用时整组自动跳过，用例文档见 TC-06-NCBI在线访问.md）"""

import os
import sys
import time
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import requests

from src.ncbi_blast import NCBIBLAST

TEST_SEQUENCE = 'AUGAGGAGAUCUCCUCAGGCA'  # 21nt，与模块内示例一致


class TestT06NCBIOnline(unittest.TestCase):
    """真实访问 NCBI BLAST 数据库"""

    @classmethod
    def setUpClass(cls):
        cls.blast = NCBIBLAST()
        # 本机直连 NCBI 稳定，但长轮询经系统代理易被中断，测试进程内绕过代理
        os.environ['NO_PROXY'] = 'blast.ncbi.nlm.nih.gov'
        os.environ['no_proxy'] = 'blast.ncbi.nlm.nih.gov'
        # 连通性探测：失败则跳过整组在线用例
        try:
            resp = requests.get(cls.blast.base_url, params={'CMD': 'Info'}, timeout=15)
            if resp.status_code != 200:
                raise unittest.SkipTest('NCBI 响应状态码 %d，跳过在线用例' % resp.status_code)
        except unittest.SkipTest:
            raise
        except Exception as e:
            raise unittest.SkipTest('无法连接 NCBI 数据库（%s），跳过在线用例' % type(e).__name__)

    def _submit(self):
        """提交真实任务并断言拿到 RID"""
        rid = self.blast.submit_blast(TEST_SEQUENCE)
        self.assertIsNotNone(rid)
        self.assertTrue(len(rid) > 0)
        return rid

    def _wait_xml(self, rid, timeout=180, poll_interval=5):
        """轮询至 COMPLETE，返回原始 XML 文本"""
        deadline = time.time() + timeout
        xml_text = None
        while time.time() < deadline:
            status, result = self.blast.check_status(rid)
            if status == 'COMPLETE' and result:
                xml_text = result
                break
            if status == 'FAILED':
                break
            time.sleep(poll_interval)
        return xml_text

    def test_tc06_001_connectivity_and_client(self):
        """TC06-001 连通性探测与客户端初始化"""
        resp = requests.get(self.blast.base_url, params={'CMD': 'Info'}, timeout=15)
        self.assertEqual(resp.status_code, 200)
        self.assertEqual(self.blast.base_url, 'https://blast.ncbi.nlm.nih.gov/Blast.cgi')

    def test_tc06_002_submit_real_task(self):
        """TC06-002 提交真实BLAST任务"""
        self._submit()

    def test_tc06_003_check_real_status(self):
        """TC06-003 查询真实任务状态"""
        rid = self._submit()
        status, _ = self.blast.check_status(rid)
        self.assertIn(status, {'WAITING', 'RUNNING', 'COMPLETE', 'FAILED'})

    def test_tc06_004_run_blast_online(self):
        """TC06-004 完整在线流程 run_blast"""
        result = self.blast.run_blast(TEST_SEQUENCE, timeout=180, poll_interval=5)
        self.assertIsNotNone(result, 'NCBI 在线完整流程未在限定时间内返回结果')
        self.assertIn('query_length', result)
        self.assertIn('hits', result)
        self.assertIsInstance(result['hits'], list)

    def test_tc06_005_parse_real_xml(self):
        """TC06-005 真实XML报文解析"""
        rid = self._submit()
        xml_text = self._wait_xml(rid)
        self.assertIsNotNone(xml_text, 'NCBI 任务未在限定时间内完成')
        result = self.blast.parse_result(xml_text)
        self.assertIsNotNone(result)
        self.assertIn('hits', result)
        # 真实报文可能没有命中（hits 为空），有命中时字段应完整
        for hit in result['hits']:
            self.assertIn('accession', hit)
            self.assertIn('identity_percent', hit)
            self.assertIn('e_value', hit)


if __name__ == '__main__':
    unittest.main(verbosity=2)
