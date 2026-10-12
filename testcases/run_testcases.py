# -*- coding: utf-8 -*-
"""测试用例执行入口：自动发现 testcases 下全部用例脚本，运行并生成 Markdown 报告

用法（在项目根目录执行）:
    python testcases/run_testcases.py
"""

import datetime
import os
import sys
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

REPORT_DIR = os.path.join(HERE, 'reports')


class CaseResultCollector(unittest.TextTestResult):
    """收集每条用例的执行结果与耗时，用于生成报告"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.records = []

    def _record(self, test, outcome, message=''):
        self.records.append({
            'case_id': test.shortDescription() or test.id(),
            'module': type(test).__module__,
            'outcome': outcome,
            'duration': self._last_duration,
            'message': message.strip(),
        })

    def _record_prefix(self):
        pass

    def startTest(self, test):
        self._start_time = time.perf_counter()
        super().startTest(test)

    def addSuccess(self, test):
        self._last_duration = time.perf_counter() - self._start_time
        super().addSuccess(test)
        self._record(test, '通过')

    def addError(self, test, err):
        self._last_duration = time.perf_counter() - self._start_time
        super().addError(test, err)
        self._record(test, '错误', self._exc_message(err))

    def addFailure(self, test, err):
        self._last_duration = time.perf_counter() - self._start_time
        super().addFailure(test, err)
        self._record(test, '失败', self._exc_message(err))

    def addSkip(self, test, reason):
        self._last_duration = time.perf_counter() - self._start_time
        super().addSkip(test, reason)
        self._record(test, '跳过', reason)

    @staticmethod
    def _exc_message(err):
        if not err:
            return ''
        return '%s: %s' % (err[0].__name__, err[1])


def build_report(records, stats, elapsed):
    """生成 Markdown 测试报告"""
    now = datetime.datetime.now()
    lines = []
    lines.append('# 测试报告')
    lines.append('')
    lines.append('- 执行时间：%s' % now.strftime('%Y-%m-%d %H:%M:%S'))
    lines.append('- Python 版本：%d.%d.%d' % sys.version_info[:3])
    lines.append('- 工作目录：%s' % ROOT)
    lines.append('- 总耗时：%.2f 秒' % elapsed)
    lines.append('- 结果汇总：共 %d 条 | 通过 %d | 失败 %d | 错误 %d | 跳过 %d' % (
        stats['total'], stats['passed'], stats['failed'], stats['errors'], stats['skipped']))
    lines.append('')

    status_mark = {'通过': 'PASS', '失败': 'FAIL', '错误': 'ERROR', '跳过': 'SKIP'}
    lines.append('| 用例编号 | 用例标题 | 所在脚本 | 状态 | 耗时(秒) | 备注 |')
    lines.append('|---------|---------|---------|------|---------|------|')
    for r in records:
        title = r['case_id']
        case_id = ''
        if title.startswith('TC'):
            parts = title.split(' ', 1)
            case_id = parts[0]
            title = parts[1] if len(parts) > 1 else ''
        note = r['message'].splitlines()[0] if r['message'] else ''
        lines.append('| %s | %s | %s | %s | %.3f | %s |' % (
            case_id, title, r['module'], status_mark[r['outcome']], r['duration'], note))
    lines.append('')

    failed = [r for r in records if r['outcome'] in ('失败', '错误')]
    if failed:
        lines.append('## 失败/错误详情')
        lines.append('')
        for r in failed:
            lines.append('### %s（%s）' % (r['case_id'], r['outcome']))
            lines.append('')
            lines.append('```')
            lines.append(r['message'])
            lines.append('```')
            lines.append('')
    return '\n'.join(lines)


def main():
    loader = unittest.TestLoader()
    suite = loader.discover(start_dir=HERE, pattern='test_tc*.py', top_level_dir=HERE)

    runner = unittest.TextTestRunner(resultclass=CaseResultCollector, verbosity=2)
    start = time.perf_counter()
    result = runner.run(suite)
    elapsed = time.perf_counter() - start

    stats = {
        'total': result.testsRun,
        'passed': len(result.records) and sum(1 for r in result.records if r['outcome'] == '通过'),
        'failed': sum(1 for r in result.records if r['outcome'] == '失败'),
        'errors': sum(1 for r in result.records if r['outcome'] == '错误'),
        'skipped': sum(1 for r in result.records if r['outcome'] == '跳过'),
    }

    os.makedirs(REPORT_DIR, exist_ok=True)
    report_name = '测试报告_%s.md' % datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    report_path = os.path.join(REPORT_DIR, report_name)
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(build_report(result.records, stats, elapsed))

    print()
    print('=' * 60)
    print('用例总数: %d | 通过: %d | 失败: %d | 错误: %d | 跳过: %d' % (
        stats['total'], stats['passed'], stats['failed'], stats['errors'], stats['skipped']))
    print('测试报告已生成: %s' % report_path)
    print('=' * 60)

    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
