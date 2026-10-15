# spec 目录说明

本目录承担两个职责：

1. 存放 PyInstaller 打包配置文件（`.spec`）
2. 作为功能变更与调试记录的跟踪中心，后续每次功能变化和 debug 都在此登记

## 目录结构

```
spec/
├── README.md                        # 本说明文件
├── PROGRESS.md                      # 开发进度总览（记录状态汇总，需与记录同步维护）
├── NucleicAcidDesigner.spec         # 打包配置（console 调试版）
├── gui.spec                         # 打包配置（windowed 发布版）
├── archive/                         # 历史 spec 文件（仅存档，不再维护）
├── changes/                         # 功能变更记录（每次功能变化一个文件）
└── debug/                           # 调试记录（每次问题排查一个文件）
```

## 打包用法

spec 文件已移入本目录，PyInstaller 以 spec 文件所在目录解析相对路径，内部路径已调整为指向项目根目录。
打包时在项目根目录执行：

```bash
# 发布版（无控制台窗口）
pyinstaller spec/gui.spec

# 调试版（带控制台，便于查看报错输出）
pyinstaller spec/NucleicAcidDesigner.spec
```

`archive/` 中的 spec 为早期命名版本，配置不完整（缺少 hiddenimports），仅作历史存档，不要再使用。

## 记录状态管理（生命周期）

每条变更/调试记录在文件头部元数据中维护一个 `状态` 字段，用统一关键词标记当前阶段，配合 `spec/PROGRESS.md` 实现进度同步管理。

### 变更记录状态（spec/changes/）

| 关键词 | 含义 | 进入条件 |
|--------|------|---------|
| plan | 已规划 | 方案已写入记录，尚未动代码 |
| developing | 开发中 | 已开始写代码，功能未完成 |
| testing | 验证中 | 功能完成，正在跑 testcases/ 或手工验证 |
| done | 已完成 | 验证通过，验证方式已登记 |
| abandoned | 已放弃 | 方案变更或需求取消，需在记录中写明原因 |

### 调试记录状态（spec/debug/）

| 关键词 | 含义 | 进入条件 |
|--------|------|---------|
| open | 排查中 | 问题确认，正在定位根因 |
| resolved | 已解决 | 根因、修复方案、验证三者齐备 |
| wontfix | 不修复 | 规避或接受现状，需写明原因 |

### 状态流转与同步规则

1. 状态沿生命周期正向流转；done / abandoned / resolved / wontfix 为终态，回退需在记录中补充说明
2. 每次状态变化必须同步更新 `spec/PROGRESS.md` 对应行，做到"记录即进度"
3. 终态记录不再修改正文，相关新内容开新文件并在开头注明关联记录
4. 任何变更进入 done 前必须已有对应的测试用例（新增或更新 testcases/），并已在"验证方式"中登记执行结果

### Agent 工作流约定（供后续 AI 会话遵循）

1. 任何开发会话开工前，先读 `spec/PROGRESS.md` 掌握当前进度与待办
2. 动代码前先在 `spec/changes/`（或 `spec/debug/`）建立记录文件，初始状态 `plan`
3. 开始编码时状态改为 `developing`，编码完成后改为 `testing`
4. 每增加新功能或修改既有功能，必须同步新增/更新 testcases/ 对应用例；无对应测试的变更不得置 done
5. 测试通过（如 `python testcases/run_testcases.py`）后改为 `done` / `resolved`，并同步更新 `spec/PROGRESS.md`
6. 代码与 spec 记录同步演进：一次功能变更对应一个记录文件与一次提交

## 功能变更跟踪

每实现一个新功能或修改既有功能，在 `changes/` 下新建记录文件：

- 命名规则：`YYYY-MM-DD-功能名.md`（例如 `2026-10-12-siRNA评分优化.md`）
- 一次功能变更一个文件，追加不改旧文件

变更记录模板（复制使用）：

```markdown
# 功能变更：<功能名>

- 日期：YYYY-MM-DD
- 状态：plan（可选值：plan / developing / testing / done / abandoned）
- 类型：新增功能 / 功能修改 / 功能删除
- 涉及模块：<如 src/sequence_design.py、gui.py>

## 变更内容

<做了什么，为什么做>

## 影响范围

<哪些模块/页面受影响，是否影响已有功能>

## 验证方式

<如何验证功能正常，例如运行了哪些命令或测试>
```

## 调试跟踪

每次排查问题（bug 修复、异常行为定位），在 `debug/` 下新建记录文件：

- 命名规则：`YYYY-MM-DD-问题简述.md`（例如 `2026-10-12-GUI导出CSV乱码.md`）
- 重点是记录根因与修复思路，避免同类问题重复排查

调试记录模板（复制使用）：

```markdown
# 调试记录：<问题简述>

- 日期：YYYY-MM-DD
- 状态：open（可选值：open / resolved / wontfix）

## 现象

<异常表现，触发步骤，报错信息>

## 根因

<定位到的根本原因>

## 修复方案

<改动了哪些文件哪些位置，如何修复>

## 验证

<复现步骤重跑结果，确认问题消失>
```
