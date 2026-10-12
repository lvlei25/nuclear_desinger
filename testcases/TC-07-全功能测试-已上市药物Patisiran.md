# TC-07 全功能测试用例：已上市药物 Patisiran（Onpattro）

以 FDA 2018 年批准的首个 siRNA 药物 Patisiran（ALN-18328，靶点 TTR mRNA 3'UTR，LNP 肝靶向，治疗 hATTR 淀粉样变性）的真实公开序列作为输入，串起工具链全部功能模块，验证端到端可用性。

公共前置（药物序列来自 IUPHAR/BPS Guide to Pharmacology 与 FDA 申报信息，去除 3' dTdT 悬垂后的 RNA 部分，各 19nt）：

- Passenger/正义链 SENSE = `AUGGAAUACUCUUGGUUAC`（真实修饰位点：pos7/17 为 2'-OMe，3' 端 dTdT）
- Guide/反义链 ANTISENSE = `GUAACCAAGAGUAUUCCAU`（真实修饰位点：pos2/12/14/15/19 为 2'-OMe，pos5/6/16/17 为 2'-OMe C，3' 端 dTdT）

| 编号 | 用例标题 | 输入 | 预期结果 |
|------|---------|------|---------|
| TC07-001 | 药物档案与序列完整性 | 两条序列字符串 | 各 19nt、仅含 AUGC、通过序列合法性校验；药物元数据（商品名/靶点/平台/年份）登记完整 |
| TC07-002 | 双链互补性验证 | reverse_complement(SENSE) | 恰好等于 ANTISENSE（真实上市药物 guide = passenger 反向互补） |
| TC07-003 | 理化性质评估 | SENSE | GC = 36.84%（处于 HF1 的 25-65% 区间）；Tm = 44.62；无免疫刺激性基序 |
| TC07-004 | 综合评估 | comprehensive_sirna_evaluation(SENSE) | final_score = 69.52、recommendation = '良好'、Reynolds = 4、Amarzguioui = 5、Ui-Tei 通过、硬过滤通过 |
| TC07-005 | LNP 修饰方案 | design_siRNA_modification(SENSE, ANTISENSE, 'LNP') | platform='LNP'、duplex='19+19'、validation.valid=True、sense 输出以 "-dTdT" 结尾（与真实药物 3' dTdT 悬垂一致）、AS 输出以 "p" 开头（5' 磷酸化） |
| TC07-006 | 递送方式与剂量 | recommend_delivery('肝脏') / calculate_dosage(19,'肝脏') | 推荐列表含 'LNP'（Patisiran 即 LNP 肝靶向）；剂量 = 0.95 mg/kg |
| TC07-007 | 全流程端到端存档 | 评估结果 + 修饰方案 + 药物档案写入 JSON/CSV 后读回 | 全链路无异常，读回字段完整（含 final_score、sense_modified、drug 元数据） |

## 行为登记（不修改源码）

- 真实药物 SENSE 在本工具综合评估中得 69.52 分（'良好'）而非 '推荐'，且热力学不对称性判据不通过：Patisiran 的实际效力依赖其真实的 2'-OMe 修饰模式与 LNP 递送，本工具的评分规则针对未修饰裸序列设计。该差异符合预期，用例按实际值断言，用于后续评分规则演进时对照。
