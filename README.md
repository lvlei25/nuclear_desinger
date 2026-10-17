# 小核酸药物设计工具包 / Nucleic Acid Drug Design Toolkit

<p align="center">
  <a href="#中文">中文</a> | <a href="#english">English</a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-1.1.0-blue" alt="version">
  <img src="https://img.shields.io/badge/python-3.9%2B-green" alt="python">
  <img src="https://img.shields.io/badge/license-MIT-lightgrey" alt="license">
</p>

---

## 中文

一个用于小核酸药物（siRNA / ASO / shRNA）设计的 Python 工具包，包含序列设计、化学修饰、数据分析与 NCBI BLAST 同源性分析等功能，并提供图形界面与命令行两种使用方式。

### 环境要求

- Python 3.9+（GUI 基于标准库 tkinter，无需额外安装）
- 依赖库见 [requirements.txt](requirements.txt)

### 安装

```bash
# 克隆仓库
git clone git@github.com:lvlei25/nuclear_desinger.git
cd nuclear_desinger

# （可选）创建虚拟环境
python -m venv .venv
# Windows
.venv\Scripts\activate
# Linux / macOS
source .venv/bin/activate

# 安装依赖
pip install -r requirements.txt
```

### 项目结构

```
nuclear_desinger/
├── src/                    # 源代码目录
│   ├── sequence_design.py      # 序列设计与评分模块
│   ├── chemical_modification.py  # 化学修饰模块
│   ├── data_processing.py      # 数据处理模块
│   ├── ncbi_blast.py           # NCBI BLAST 同源性分析模块
│   ├── logger.py               # 日志模块（按天落盘 log/ 目录）
│   └── __init__.py             # 包初始化文件
├── data/                   # 数据目录
│   └── sequences/              # 序列数据库
├── config/                 # 配置目录
│   └── settings.py             # 配置文件（含 ScoringConfig 评分参数）
├── spec/                   # 变更与调试记录（spec 跟踪机制）
├── testcases/              # 测试用例体系（9 组 98 条，含报告生成）
├── tests/                  # 早期单元测试（13 条）
├── log/                    # 运行日志输出目录（按天生成）
├── gui.py                  # 图形用户界面
├── demo.py                 # 命令行演示脚本
├── CHANGELOG.md            # 版本历史
└── README.md               # 项目说明
```

### 使用方式

#### 方式一：图形用户界面（推荐）

```bash
python gui.py
```

GUI 功能特点：

- **序列设计页面**：输入 mRNA 序列，设计 siRNA，计算 GC 含量和 Tm 值，查找基序
- **化学修饰页面**：添加和管理化学修饰，获取推荐修饰策略
- **递送系统页面**：根据目标组织推荐递送方式和剂量
- **序列数据库页面**：保存、加载和导出序列数据

#### 方式二：命令行演示

```bash
python demo.py
```

### 功能模块

#### 1. 序列设计 (sequence_design.py)

- 互补链和反向互补链计算
- GC 含量计算
- siRNA 序列设计与综合评估（Reynolds / Ui-Tei / Amarzguioui / 热力学不对称性四套评分规则 + 硬性过滤）
- ASO（反义寡核苷酸）设计
- shRNA 序列生成
- 脱靶效应检测与特异性评分
- 评分权重与阈值支持自定义配置（`config/settings.py` 的 `ScoringConfig`）

#### 2. 化学修饰 (chemical_modification.py)

- 多种修饰类型支持（2OMe, PS, LNA, FANA 等）
- 修饰策略推荐
- 递送系统设计
- 剂量计算

#### 3. 数据处理 (data_processing.py)

- FASTA 文件读写
- CSV 文件处理
- JSON 数据处理
- 序列数据库管理

#### 4. 同源性分析 (ncbi_blast.py)

- 基于 NCBI BLAST API 的序列同源性检索
- 异步提交与结果轮询
- XML 结果解析与脱靶位点评估

### 测试运行

```bash
# 全量测试用例（自动生成 Markdown 报告到 testcases/reports/）
python testcases/run_testcases.py

# 跳过 NCBI 在线用例（CI / 无网环境）
# PowerShell:  $env:SKIP_ONLINE_TESTS='1'; python testcases/run_testcases.py
# Bash:        SKIP_ONLINE_TESTS=1 python testcases/run_testcases.py

# 早期单元测试
python -m pytest tests/ -v
```

推送 / PR 时 GitHub Actions 自动运行全量用例（见 [.github/workflows/tests.yml](.github/workflows/tests.yml)）。

### 日志

程序运行日志按天写入 `log/app_YYYYMMDD.log`（UTF-8），GUI 未捕获异常也会落盘，便于无人值守场景事后排查。

### 许可证

MIT License

### 作者

小核酸药物设计团队 / Nucleic Acid Drug Design Team

---

## English

A Python toolkit for small nucleic acid drug design (siRNA / ASO / shRNA), covering sequence design, chemical modification, data processing and NCBI BLAST homology analysis, with both GUI and command-line interfaces.

### Requirements

- Python 3.9+ (the GUI is based on the standard-library tkinter, no extra installation needed)
- Dependencies: see [requirements.txt](requirements.txt)

### Installation

```bash
# Clone the repository
git clone git@github.com:lvlei25/nuclear_desinger.git
cd nuclear_desinger

# (Optional) Create a virtual environment
python -m venv .venv
# Windows
.venv\Scripts\activate
# Linux / macOS
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### Project Structure

```
nuclear_desinger/
├── src/                    # Source code
│   ├── sequence_design.py      # Sequence design & scoring
│   ├── chemical_modification.py  # Chemical modification
│   ├── data_processing.py      # Data processing
│   ├── ncbi_blast.py           # NCBI BLAST homology analysis
│   ├── logger.py               # Logging (daily files under log/)
│   └── __init__.py             # Package init
├── data/                   # Data directory
│   └── sequences/              # Sequence database
├── config/                 # Configuration
│   └── settings.py             # Settings (incl. ScoringConfig scoring parameters)
├── spec/                   # Change & debug records (spec tracking)
├── testcases/              # Test suite (9 groups / 98 cases, with report generation)
├── tests/                  # Early unit tests (13 cases)
├── log/                    # Runtime log output (generated daily)
├── gui.py                  # Graphical user interface
├── demo.py                 # Command-line demo script
├── CHANGELOG.md            # Version history
└── README.md               # This file
```

### Usage

#### Option 1: Graphical User Interface (recommended)

```bash
python gui.py
```

GUI features:

- **Sequence Design**: input an mRNA sequence, design siRNA, compute GC content and Tm, find motifs
- **Chemical Modification**: add and manage modifications, get recommended modification strategies
- **Delivery System**: get delivery method and dosage recommendations by target tissue
- **Sequence Database**: save, load and export sequence data

#### Option 2: Command-line Demo

```bash
python demo.py
```

### Feature Modules

#### 1. Sequence Design (sequence_design.py)

- Complement and reverse-complement computation
- GC content calculation
- siRNA design with comprehensive evaluation (Reynolds / Ui-Tei / Amarzguioui / thermodynamic asymmetry scoring rules + hard filters)
- ASO (antisense oligonucleotide) design
- shRNA sequence generation
- Off-target detection and specificity scoring
- Customizable scoring weights and thresholds (`ScoringConfig` in `config/settings.py`)

#### 2. Chemical Modification (chemical_modification.py)

- Multiple modification types (2OMe, PS, LNA, FANA, etc.)
- Modification strategy recommendation
- Delivery system design
- Dosage calculation

#### 3. Data Processing (data_processing.py)

- FASTA file I/O
- CSV file handling
- JSON data processing
- Sequence database management

#### 4. Homology Analysis (ncbi_blast.py)

- Sequence homology search via the NCBI BLAST API
- Asynchronous submission and result polling
- XML result parsing and off-target site evaluation

### Running Tests

```bash
# Full test suite (auto-generates a Markdown report under testcases/reports/)
python testcases/run_testcases.py

# Skip NCBI online cases (CI / offline environments)
SKIP_ONLINE_TESTS=1 python testcases/run_testcases.py

# Early unit tests
python -m pytest tests/ -v
```

GitHub Actions runs the full suite on push / pull requests (see [.github/workflows/tests.yml](.github/workflows/tests.yml)).

### Logging

Runtime logs are written daily to `log/app_YYYYMMDD.log` (UTF-8). Uncaught GUI exceptions are also logged, enabling post-mortem debugging in unattended scenarios.

### License

MIT License

### Authors

Nucleic Acid Drug Design Team
