# 小核酸药物设计工具包

一个用于小核酸药物设计的Python工具包，包含序列设计、化学修饰、数据分析等功能。

## 环境要求

- Python 3.9+（GUI 基于标准库 tkinter，无需额外安装）
- 依赖库见 [requirements.txt](requirements.txt)

## 安装

```bash
# 克隆仓库
git clone https://github.com/<your-username>/<your-repo>.git
cd 小核酸药物设计

# （可选）创建虚拟环境
python -m venv .venv
# Windows
.venv\Scripts\activate
# Linux / macOS
source .venv/bin/activate

# 安装依赖
pip install -r requirements.txt
```

## 项目结构

```
小核酸药物设计/
├── src/                    # 源代码目录
│   ├── sequence_design.py    # 序列设计模块
│   ├── chemical_modification.py  # 化学修饰模块
│   ├── data_processing.py     # 数据处理模块
│   ├── ncbi_blast.py          # NCBI BLAST 同源性分析模块
│   └── __init__.py            # 包初始化文件
├── data/                   # 数据目录
│   └── sequences/             # 序列数据库
├── config/                 # 配置目录
│   └── settings.py            # 配置文件
├── tests/                  # 测试目录
│   ├── test_sequence_design.py      # 序列设计测试
│   └── test_chemical_modification.py  # 化学修饰测试
├── docs/                   # 文档目录
├── gui.py                  # 图形用户界面（交互对话框）
├── demo.py                 # 命令行演示脚本
└── README.md               # 项目说明
```

## 使用方式

### 方式一：图形用户界面（推荐）

直接运行GUI应用程序，享受交互式对话框体验：

```bash
cd 小核酸药物设计
python gui.py
```

GUI功能特点：
- **序列设计页面**：输入mRNA序列，设计siRNA，计算GC含量和Tm值，查找基序
- **化学修饰页面**：添加和管理化学修饰，获取推荐修饰策略
- **递送系统页面**：根据目标组织推荐递送方式和剂量
- **序列数据库页面**：保存、加载和导出序列数据

### 方式二：命令行演示

```bash
cd 小核酸药物设计
python demo.py
```

## 功能模块

### 1. 序列设计 (sequence_design.py)
- 互补链和反向互补链计算
- GC含量计算
- siRNA序列设计
- ASO（反义寡核苷酸）设计
- shRNA序列生成
- 脱靶效应检测

### 2. 化学修饰 (chemical_modification.py)
- 多种修饰类型支持（2OMe, PS, LNA, FANA等）
- 修饰策略推荐
- 递送系统设计
- 剂量计算

### 3. 数据处理 (data_processing.py)
- FASTA文件读写
- CSV文件处理
- JSON数据处理
- 序列数据库管理

### 4. 同源性分析 (ncbi_blast.py)
- 基于 NCBI BLAST API 的序列同源性检索
- 异步提交与结果轮询
- XML 结果解析与脱靶位点评估

## GUI界面预览

图形用户界面包含4个主要标签页：

1. **序列设计**：输入mRNA序列，设置设计参数，一键生成siRNA
2. **化学修饰**：选择修饰类型和位置，查看推荐策略
3. **递送系统**：选择目标组织，获取递送方式和剂量推荐
4. **序列数据库**：管理保存的序列，导出为CSV格式

## 测试运行

```bash
python -m pytest tests/ -v
```

## 许可证

MIT License

## 作者

小核酸药物设计团队