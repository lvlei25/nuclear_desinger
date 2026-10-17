# CHANGELOG

本文件记录小核酸药物设计工具包的版本历史。版本号定义于 `config/settings.py` 与 `src/__init__.py` 的 `__version__`（两处保持同步）；
每次发版的详细变更见 `spec/changes/` 对应记录，形成"版本 → 变更记录 → 测试报告"追溯链。

## [1.1.0] - 2026-10-17

### 新增
- P1 工程质量：配置外置——评分权重与阈值集中到 `config/settings.py` 的 `ScoringConfig`（默认值与历史硬编码一致），`NucleicAcidDesigner` 支持注入自定义配置（spec/changes/2026-10-17-P1工程质量-配置外置-版本号-CI.md）
- CI：GitHub Actions 工作流 `.github/workflows/tests.yml`，push/PR 自动跑 testcases/ 全量用例（`SKIP_ONLINE_TESTS=1` 跳过 NCBI 在线组）
- 测试组 TC-09 配置外置（9 条用例）；TC-06 在线用例支持 `SKIP_ONLINE_TESTS` 环境变量强制跳过
- 版本追溯链：`__version__` + 本 CHANGELOG

### 修复
- GUI 基因变量未定义缺陷（spec/debug/2026-10-15-GUI基因变量未定义.md）

## [1.0.0] - 2026-10-15

### 新增
- log 日志系统：`src/logger.py` 按天落盘 `log/` 目录，入口与错误路径接入（spec/changes/2026-10-15-log日志系统.md）
- 基因别名自动解析（spec/changes/2026-10-15-基因别名自动解析.md）
- 测试组 TC-08 日志系统

### 修复
- P0 问题修复（spec/changes/2026-10-15-P0问题修复.md）
- TC06 在线用例提交重试与超时提升（spec/changes/2026-10-15-TC06在线用例提交重试.md）

## [0.9.0] - 2026-10-12

### 新增
- 初始功能包：序列设计（siRNA/ASO/shRNA，含 Reynolds/Ui-Tei/Amarzguioui/热力学不对称性四套评分）、化学修饰与递送设计、数据读写与序列数据库、NCBI BLAST 客户端、Tkinter GUI（spec/changes/2026-10-12-功能基线.md）
- 测试用例体系 testcases/（8 组）与全功能回归基准（已上市药物 Patisiran）
- spec 机制：spec/changes 变更记录 + spec/debug 调试记录 + PROGRESS.md 进度总览
