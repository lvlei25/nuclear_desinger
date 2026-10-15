# -*- coding: utf-8 -*-
"""
日志系统模块

统一配置项目日志：按天写入 log/app_YYYYMMDD.log（UTF-8），同时保留控制台输出，
便于人工查看程序运行状态与事后排查。
"""

import logging
import os
import sys
from datetime import datetime

# 日志目录（项目根目录下的 log 文件夹）
LOG_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'log')

# 已初始化标记，保证 setup_logging 幂等
_configured = False


def setup_logging(level=logging.INFO):
    """配置项目日志：文件（按天）+ 控制台双通道，重复调用幂等

    参数:
        level: 日志级别，默认 INFO
    """
    global _configured
    if _configured:
        return

    # 确保日志目录存在
    os.makedirs(LOG_DIR, exist_ok=True)

    # 日志格式：时间 [级别] 模块名: 消息
    fmt = logging.Formatter('%(asctime)s [%(levelname)s] %(name)s: %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')

    # 文件通道：按天一个文件，UTF-8 编码
    log_file = os.path.join(LOG_DIR, 'app_%s.log' % datetime.now().strftime('%Y%m%d'))
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setFormatter(fmt)
    file_handler.setLevel(level)

    # 控制台通道：保持原有控制台可见性
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setFormatter(fmt)
    console_handler.setLevel(level)

    # 挂载到根 logger，项目内所有子 logger 均可输出
    root = logging.getLogger()
    root.setLevel(level)
    root.addHandler(file_handler)
    root.addHandler(console_handler)

    _configured = True

    logging.getLogger(__name__).info('日志系统初始化完成，日志文件: %s', log_file)


def get_logger(name):
    """获取指定名称的子 logger（setup_logging 之前调用也能安全创建）"""
    return logging.getLogger(name)


def install_excepthook():
    """注册全局异常钩子：未捕获异常写入日志（GUI 无人值守时可事后排查）"""
    original_hook = sys.excepthook

    def _hook(exc_type, exc_value, exc_tb):
        if issubclass(exc_type, KeyboardInterrupt):
            original_hook(exc_type, exc_value, exc_tb)
            return
        logging.getLogger(__name__).error('未捕获异常: %s', exc_value,
                                          exc_info=(exc_type, exc_value, exc_tb))
        original_hook(exc_type, exc_value, exc_tb)

    sys.excepthook = _hook
