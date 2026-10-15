# -*- coding: utf-8 -*-
"""TC-04 同源性分析模块测试用例（requests 全部 Mock，不联网，用例文档见 TC-04-同源性分析BLAST.md）"""

import os
import sys
import unittest
from unittest.mock import MagicMock, patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.ncbi_blast import NCBIBLAST

# 构造的最小 BLAST XML（覆盖解析器实际读取的全部标签；不含 XML 声明，便于与状态标记拼接）
SAMPLE_XML = """<BlastOutput>
  <Query-def>Query_1</Query-def>
  <Query-len>21</Query-len>
  <Hit>
    <Hit_num>1</Hit_num>
    <Hit_id>gi|999|ref|NM_000001.1|</Hit_id>
    <Hit_def>synthetic construct targetgene mRNA</Hit_def>
    <Hit_accession>NM_000001</Hit_accession>
    <Hit_len>21</Hit_len>
    <Hsp>
      <Hsp_identity>19</Hsp_identity>
      <Hsp_align-len>21</Hsp_align-len>
      <Hsp_evalue>0.001</Hsp_evalue>
      <Hsp_score>39.5</Hsp_score>
      <Hsp_strand>Plus/Plus</Hsp_strand>
    </Hsp>
  </Hit>
</BlastOutput>"""

# 真实 NCBI 报文形态的 XML（BlastOutput_iterations 包裹，查询标签为 Iteration_query-def / Iteration_query-len）
REALFORM_XML = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-def>Query_1</Iteration_query-def>
      <Iteration_query-len>21</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gi|999|ref|NM_000001.1|</Hit_id>
          <Hit_def>synthetic construct targetgene mRNA</Hit_def>
          <Hit_accession>NM_000001</Hit_accession>
          <Hit_len>21</Hit_len>
          <Hsp>
            <Hsp_identity>19</Hsp_identity>
            <Hsp_align-len>21</Hsp_align-len>
            <Hsp_evalue>0.001</Hsp_evalue>
            <Hsp_score>39.5</Hsp_score>
            <Hsp_strand>Plus/Plus</Hsp_strand>
          </Hsp>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""


def make_response(status_code=200, text=''):
    """构造 Mock HTTP 响应"""
    resp = MagicMock()
    resp.status_code = status_code
    resp.text = text
    return resp


class TestT04SubmitAndPoll(unittest.TestCase):
    """提交与轮询"""

    def setUp(self):
        self.blast = NCBIBLAST()

    @patch('src.ncbi_blast.requests.post')
    def test_tc04_001_submit_success(self, mock_post):
        """TC04-001 提交成功提取RID"""
        mock_post.return_value = make_response(200, 'RID = ABC123\nRTOE = 10\n')
        self.assertEqual(self.blast.submit_blast('AUGCAUGC'), 'ABC123')

    @patch('src.ncbi_blast.requests.post')
    def test_tc04_002_submit_bad_status(self, mock_post):
        """TC04-002 提交失败状态码"""
        mock_post.return_value = make_response(500, 'error')
        self.assertIsNone(self.blast.submit_blast('AUGCAUGC'))

    @patch('src.ncbi_blast.requests.post', side_effect=ConnectionError('网络不可用'))
    def test_tc04_003_submit_exception(self, mock_post):
        """TC04-003 提交网络异常返回None"""
        self.assertIsNone(self.blast.submit_blast('AUGCAUGC'))

    @patch('src.ncbi_blast.requests.get')
    def test_tc04_004_check_status_branches(self, mock_get):
        """TC04-004 轮询状态分支"""
        cases = {
            'xxx Status=WAITING xxx': ('WAITING', None),
            'xxx Status=RUNNING xxx': ('RUNNING', None),
            'xxx Status=COMPLETE xxx': ('COMPLETE', 'xxx Status=COMPLETE xxx'),
            'no status here': ('FAILED', None),
        }
        for text, expected in cases.items():
            mock_get.return_value = make_response(200, text)
            self.assertEqual(self.blast.check_status('RID1'), expected)


class TestT04ParseAndRisk(unittest.TestCase):
    """XML 解析与脱靶风险评估"""

    def setUp(self):
        self.blast = NCBIBLAST()

    def test_tc04_005_parse_result(self):
        """TC04-005 XML结果解析"""
        result = self.blast.parse_result(SAMPLE_XML)
        self.assertIsNotNone(result)
        self.assertEqual(result['query_id'], 'Query_1')
        self.assertEqual(result['query_length'], 21)
        self.assertEqual(len(result['hits']), 1)
        hit = result['hits'][0]
        self.assertEqual(hit['accession'], 'NM_000001')
        self.assertAlmostEqual(hit['identity_percent'], 90.476, places=3)
        self.assertAlmostEqual(hit['e_value'], 0.001, places=6)
        self.assertEqual(hit['strand'], 'Plus/Plus')

    def test_tc04_006_parse_invalid_xml(self):
        """TC04-006 非法XML解析返回None"""
        self.assertIsNone(self.blast.parse_result('not xml'))

    def _risk_with_hits(self, hits, target_gene=''):
        return self.blast.get_off_target_risk({'hits': hits}, target_gene=target_gene)

    def test_tc04_007_off_target_levels(self):
        """TC04-007 脱靶风险三级分级"""
        hits = [
            {'accession': 'A', 'description': 'gene A', 'identity_percent': 95.2, 'e_value': 0.001},
            {'accession': 'B', 'description': 'gene B', 'identity_percent': 86.0, 'e_value': 0.01},
            {'accession': 'C', 'description': 'gene C', 'identity_percent': 81.0, 'e_value': 0.1},
        ]
        risk = self._risk_with_hits(hits)
        self.assertTrue(risk['has_off_target'])
        self.assertEqual(len(risk['high_risk']), 1)
        self.assertEqual(len(risk['medium_risk']), 1)
        self.assertEqual(len(risk['low_risk']), 1)
        self.assertEqual(risk['total_hits'], 3)

    def test_tc04_008_below_threshold_ignored(self):
        """TC04-008 低于阈值过滤"""
        hits = [
            {'accession': 'A', 'description': 'gene A', 'identity_percent': 95.2, 'e_value': 0.001},
            {'accession': 'D', 'description': 'gene D', 'identity_percent': 70.0, 'e_value': 1.0},
        ]
        risk = self._risk_with_hits(hits)
        self.assertEqual(risk['total_hits'], 1)
        self.assertEqual(len(risk['high_risk']), 1)

    def test_tc04_009_target_gene_excluded(self):
        """TC04-009 目标基因自身排除"""
        hits = [
            {'accession': 'E', 'description': 'Homo sapiens targetgene mRNA', 'identity_percent': 95.0, 'e_value': 0.001},
        ]
        risk = self._risk_with_hits(hits, target_gene='targetgene')
        self.assertFalse(risk['has_off_target'])
        self.assertTrue(risk['target_excluded'])
        self.assertEqual(risk['total_hits'], 0)

    @patch('src.ncbi_blast.requests.get')
    @patch('src.ncbi_blast.requests.post')
    def test_tc04_010_run_blast_full_flow(self, mock_post, mock_get):
        """TC04-010 run_blast完整流程（首轮即完成）"""
        mock_post.return_value = make_response(200, 'RID = ABC123\n')
        # 状态标记置于 XML 注释内，保证整段响应文本可被 XML 解析
        mock_get.return_value = make_response(200, '<?xml version="1.0"?>\n<!-- Status=COMPLETE -->\n' + SAMPLE_XML)
        result = self.blast.run_blast('AUGCAUGC')
        self.assertIsNotNone(result)
        self.assertEqual(result['query_id'], 'Query_1')
        self.assertEqual(len(result['hits']), 1)

    def test_tc04_011_realform_xml_parse(self):
        """TC04-011 真实 NCBI 形态 XML 解析（P0-3：Iteration_query-def/len 标签兼容）"""
        result = self.blast.parse_result(REALFORM_XML)
        self.assertIsNotNone(result)
        self.assertEqual(result['query_id'], 'Query_1')
        self.assertEqual(result['query_length'], 21)
        self.assertEqual(len(result['hits']), 1)
        self.assertEqual(result['hits'][0]['accession'], 'NM_000001')


if __name__ == '__main__':
    unittest.main(verbosity=2)
