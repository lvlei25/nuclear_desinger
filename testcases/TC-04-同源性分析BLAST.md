# TC-04 同源性分析模块用例（src/ncbi_blast.py）

公共前置：全部用例通过 `unittest.mock.patch` 替换 `src.ncbi_blast.requests` 的 post/get 响应，不发起真实网络请求。

| 编号 | 用例标题 | 输入 | 预期结果 |
|------|---------|------|---------|
| TC04-001 | 提交成功提取 RID | post 200，响应含 "RID = ABC123" | 返回 "ABC123" |
| TC04-002 | 提交失败状态码 | post 返回 500 | 返回 None |
| TC04-003 | 提交网络异常 | post 抛 ConnectionError | 捕获异常返回 None |
| TC04-004 | 轮询状态分支 | get 响应分别含 Status=WAITING/RUNNING/COMPLETE/其他 | (WAITING,None) / (RUNNING,None) / (COMPLETE,响应文本) / (FAILED,None) |
| TC04-005 | XML 结果解析 | 构造含 1 个 Hit（19/21）的 XML | query_id="Query_1"、query_length=21、hits[0].identity_percent 约 90.48、e_value/score/strand 正确 |
| TC04-006 | 非法 XML 解析 | "not xml" | 捕获异常返回 None |
| TC04-007 | 脱靶风险三级分级 | identity 95.2/86.0/81.0 三条命中 | 分别进入 high_risk/medium_risk/low_risk，total_hits=3，has_off_target=True |
| TC04-008 | 低于阈值过滤 | 追加 identity 70.0 命中 | 不计入任何风险级别，total_hits 不变 |
| TC04-009 | 目标基因自身排除 | identity 95 且描述含 "targetgene" | 不计入风险，target_excluded=True |
| TC04-010 | run_blast 完整流程 | Mock 提交成功且首轮即 COMPLETE | 返回解析后的结果字典（含 hits），全程无真实网络与长轮询 |
| TC04-011 | 真实形态 XML 解析 | parse_result(REALFORM_XML)（Iteration_query-def/len 标签） | query_id="Query_1"、query_length=21、hits 正常解析（P0-3 兼容修复） |
