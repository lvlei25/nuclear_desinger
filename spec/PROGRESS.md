# 开发进度总览（PROGRESS）

> 本文件是 spec 记录状态的唯一汇总视图：任何记录的状态变化必须同步更新本表，做到"记录即进度"。
> 状态关键词定义见 `spec/README.md` 的"记录状态管理（生命周期）"。
> 最后更新：2026-10-17

## 模块状态总览

| 模块 | 说明 | 状态 |
|------|------|------|
| src/sequence_design.py | 序列设计与评分（siRNA/ASO/shRNA） | 稳定（P0-4 边界统一：空/短序列抛带明确信息的 ValueError） |
| src/chemical_modification.py | 化学修饰与递送设计 | 稳定 |
| src/data_processing.py | 数据读写与序列数据库 | 稳定 |
| src/ncbi_blast.py | NCBI BLAST 客户端 | 稳定（真实访问缺陷已修复；parse_result 已兼容真实 XML 标签） |
| gui.py | Tkinter 图形界面 | 稳定（错误路径已接入 logging） |
| src/logger.py | 日志系统（按天落盘 log/ 目录） | 稳定（2026-10-15 新增） |
| config/settings.py | ScoringConfig 评分配置（外置权重与阈值） | 稳定（2026-10-17 外置，默认值=历史硬编码） |
| testcases/ | 测试用例体系（9 组 98 条） | 持续维护 |
| tests/ | 早期单元测试（13 条） | 全绿（2026-10-15 修正 4 处历史断言） |

## 变更记录状态（spec/changes/）

| 记录 | 状态 | 完成日期 |
|------|------|---------|
| 2026-10-12-功能基线.md | done | 2026-10-12 |
| 2026-10-12-测试用例体系.md | done | 2026-10-12 |
| 2026-10-12-全功能测试-已上市药物Patisiran.md | done | 2026-10-12 |
| 2026-10-15-P0问题修复.md | done | 2026-10-15 |
| 2026-10-15-log日志系统.md | done | 2026-10-15 |
| 2026-10-15-基因别名自动解析.md | done | 2026-10-15 |
| 2026-10-15-TC06在线用例提交重试.md | done | 2026-10-15 |
| 2026-10-17-P1工程质量-配置外置-版本号-CI.md | done | 2026-10-17 |

## 调试记录状态（spec/debug/）

| 记录 | 状态 | 解决日期 |
|------|------|---------|
| 2026-10-12-NCBI任务完成状态判断.md | resolved | 2026-10-12 |
| 2026-10-15-GUI基因变量未定义.md | resolved | 2026-10-15 |

## 待办池（2026-10-15 项目评审归档，开工时先建 plan 记录）

### P1 工程质量

- ~~配置外置~~：已完成（2026-10-17，见 spec/changes/2026-10-17-P1工程质量-配置外置-版本号-CI.md，ScoringConfig 外置全部权重与阈值）
- ~~版本号~~：已完成（2026-10-17，__version__=1.1.0 + CHANGELOG.md 追溯链）
- ~~CI~~：已完成（2026-10-17，.github/workflows/tests.yml，SKIP_ONLINE_TESTS=1 跳过 TC-06 在线组）
- GUI 业务逻辑抽离：gui.py 2000+ 行单文件，界面与业务耦合，抽成 src/ 可测模块
- logging 剩余范围：基础版已完成（log/ 目录按天落盘、入口与错误路径接入）；剩余为剩余散落输出的清理与 GUI 界面内嵌日志查看（可选）
- 已上市药物知识库：扩展 data/approved_drugs.json（Inclisiran/Givosiran/Lumasiran/Vutrisiran 等），作回归基准与修饰参考

### P2 科学性增强

- 种子区特异性评估：第 2-8 位种子区 vs 3'UTR 匹配评分，配合本地 BLAST 库（makeblastdb）摆脱联网依赖
- 热力学 profile：全序列 ΔG 分布曲线（不对称性 + 终端稳定性），替代仅前 4 位判据
- 修饰效果量化：2'-OMe/F/PS 对免疫/稳定性/脱靶的修正系数，使"修饰后评分"成为可能（Patisiran 用例为对照基准）
- miRNA/ASO 完整设计链路：目前 suggest_modifications 支持 miRNA 但无完整设计流程
