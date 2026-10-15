# testcases 测试用例说明

本目录是功能测试用例的统一管理入口，与 `tests/`（开发者单元测试）相互独立。

## 目录结构

```
testcases/
├── README.md                        # 本说明
├── TC-01-序列设计.md                # 用例文档（预期结果）
├── TC-02-化学修饰.md
├── TC-03-数据处理.md
├── TC-04-同源性分析BLAST.md
├── TC-05-综合评估与GUI.md
├── TC-06-NCBI在线访问.md
├── TC-07-全功能测试-已上市药物Patisiran.md
├── TC-08-日志系统.md
├── test_tc01_sequence_design.py     # 用例对应的可执行脚本
├── test_tc02_chemical_modification.py
├── test_tc03_data_processing.py
├── test_tc04_ncbi_blast.py
├── test_tc05_integration_gui.py
├── test_tc06_ncbi_online.py
├── test_tc07_marketed_drug.py
├── test_tc08_logging.py
├── run_testcases.py                 # 执行入口（自动发现用例并生成报告）
└── reports/                         # 测试报告输出目录（运行后生成）
```

## 测试范围

| 编号 | 模块 | 覆盖内容 | 用例数 |
|------|------|---------|--------|
| TC-01 | src/sequence_design.py | 链计算、GC/Tm、siRNA 设计、四套评分规则、硬过滤、脱靶、候选查找、ASO/shRNA/优化 | 31 |
| TC-02 | src/chemical_modification.py | 修饰增删改查、标记/合成格式、修饰方案设计（GalNAc/LNP）、递送与剂量 | 14 |
| TC-03 | src/data_processing.py | FASTA/CSV/JSON 读写、目标基因加载、序列数据库 | 8 |
| TC-04 | src/ncbi_blast.py | 提交/轮询/解析/脱靶风险评估（requests 全部 Mock，不联网） | 11 |
| TC-05 | 综合评估 + gui.py + demo.py | 综合评分、端到端流程、demo 冒烟、GUI 导入冒烟、基因别名解析（Mock） | 9 |
| TC-06 | src/ncbi_blast.py（在线集成） | 真实访问 NCBI 数据库：提交/状态/完整流程/真实 XML 解析/别名解析；网络不可用时整组自动跳过 | 6 |
| TC-07 | 全功能（已上市药物 Patisiran） | 用 FDA 获批药物 Onpattro 的真实双链序列串起设计/评估/修饰/递送/存档全流程 | 7 |
| TC-08 | src/logger.py + 错误路径接入 | 日志初始化幂等、按天文件落盘（UTF-8）、格式规范、ncbi_blast 错误落盘 | 5 |

合计 91 条用例，按"功能点 + 边界值 + 异常输入"三层设计；TC-06 依赖外网，离线环境仅影响该组（跳过）。

## 运行环境

- Python 3.9+（验证环境为 3.13），仅依赖标准库 unittest + 项目已有依赖（requests）
- BLAST 模块通过 unittest.mock 模拟网络响应，运行全程无需外网

## 执行方式

在项目根目录执行：

```bash
# 全部用例并生成报告
python testcases/run_testcases.py

# 单独运行某一组
python -m unittest testcases.test_tc01_sequence_design -v
```

运行结束后在 `reports/` 下生成 Markdown 测试报告（含每条用例的通过状态与耗时）。

## 用例文档约定

- 每条用例编号 `TCxx-yyy`，与测试脚本中 unittest 方法的 docstring 一一对应
- 文档记录"前置条件、输入、预期结果"；实际执行结果统一记录在 `reports/` 报告中
- 边界实现行为按登记演进：空序列 GC 计算已由 ZeroDivisionError 统一为 ValueError（2026-10-15 P0 修复），后续行为调整均在 spec 变更记录中评估登记

## 与既有 tests/ 的关系

`tests/` 是早期单元测试，其中 4 处断言与当前实现不符（详见 TC-01 文档备注），已于 2026-10-15 的 P0 修复中修正测试数据，恢复全绿基线（实现行为未改动）。
