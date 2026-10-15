# -*- coding: utf-8 -*-
"""TC-08 日志系统测试用例（用例文档见 TC-08-日志系统.md）"""

import logging
import os
import sys
import unittest
from datetime import datetime
from unittest.mock import patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.logger import setup_logging, get_logger, install_excepthook, LOG_DIR


def today_log_file():
    """获取今天的日志文件路径"""
    return os.path.join(LOG_DIR, 'app_%s.log' % datetime.now().strftime('%Y%m%d'))


def read_log_text():
    """读取今天的日志全文"""
    with open(today_log_file(), 'r', encoding='utf-8') as f:
        return f.read()


class TestT08LoggingSetup(unittest.TestCase):
    """日志初始化与文件落盘"""

    def test_tc08_001_setup_idempotent(self):
        """TC08-001 setup_logging 幂等：重复调用不重复挂载 handler"""
        setup_logging()
        setup_logging()
        file_handlers = [h for h in logging.getLogger().handlers
                         if isinstance(h, logging.FileHandler)]
        self.assertEqual(len(file_handlers), 1)

    def test_tc08_002_log_file_created(self):
        """TC08-002 日志文件按天创建于 log/ 目录"""
        setup_logging()
        self.assertTrue(os.path.isdir(LOG_DIR))
        self.assertTrue(os.path.isfile(today_log_file()))
        self.assertIn('app_%s.log' % datetime.now().strftime('%Y%m%d'), today_log_file())

    def test_tc08_003_utf8_chinese_content(self):
        """TC08-003 中文消息以 UTF-8 写入且不出现乱码"""
        logger = get_logger('tc08.charset')
        marker = '中文写入测试TC08_运行状态正常'
        logger.info(marker)
        text = read_log_text()
        self.assertIn(marker, text)


class TestT08LoggingContent(unittest.TestCase):
    """日志内容格式与模块接入"""

    def test_tc08_004_format_and_module_name(self):
        """TC08-004 日志行含时间戳、级别与子模块名"""
        setup_logging()
        logger = get_logger('tc08.format')
        marker = '格式校验TC08_004'
        logger.warning(marker)
        text = read_log_text()
        # 同一天多次运行用例时日志会累计，匹配数 >=1 即可，格式校验取最新一条
        matched = [ln for ln in text.splitlines() if marker in ln]
        self.assertGreaterEqual(len(matched), 1)
        line = matched[-1]
        # 行首为日期时间（YYYY-MM-DD HH:MM:SS）
        self.assertRegex(line, r'^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}')
        self.assertIn('[WARNING]', line)
        self.assertIn('tc08.format', line)

    def test_tc08_005_ncbi_error_logged_to_file(self):
        """TC08-005 ncbi_blast 网络异常经 logger 落盘"""
        setup_logging()
        from src.ncbi_blast import NCBIBLAST
        blast = NCBIBLAST()
        with patch('src.ncbi_blast.requests.post', side_effect=Exception('TC08模拟网络中断')):
            rid = blast.submit_blast('AUGCAUGC')
        self.assertIsNone(rid)
        text = read_log_text()
        self.assertIn('BLAST提交失败: TC08模拟网络中断', text)
        self.assertIn('src.ncbi_blast', text)


if __name__ == '__main__':
    unittest.main(verbosity=2)
