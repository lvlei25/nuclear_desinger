# TC-06 NCBI 在线访问集成用例（src/ncbi_blast.py 真实联网）

与 TC-04（纯 Mock）互补，本组用例真实访问 NCBI BLAST 数据库，验证端到端可用性与真实报文兼容性。

公共前置：
- 本机能直连 `https://blast.ncbi.nlm.nih.gov`（实测 200，约 2 秒）
- `setUpClass` 先做连通性探测，探测失败（无网络/被防火墙拦截）时整组自动跳过（SKIP），不影响其余用例
- 在线用例耗时受 NCBI 排队影响，单条完整流程预留 180 秒超时

| 编号 | 用例标题 | 输入 | 预期结果 |
|------|---------|------|---------|
| TC06-001 | 连通性探测与客户端初始化 | GET CMD=Info | 状态码 200，NCBIBLAST 实例 base_url 正确 |
| TC06-002 | 提交真实 BLAST 任务 | 21nt 测试序列 AUGAGGAGAUCUCCUCAGGCA | 返回非空 RID 字符串 |
| TC06-003 | 查询真实任务状态 | 已获取的 RID | 状态 ∈ {WAITING, RUNNING, COMPLETE, FAILED} |
| TC06-004 | 完整在线流程 run_blast | 同上序列，timeout=600，轮询间隔 5s（NCBI 高峰期任务排队可超过 5 分钟，2026-10-15 由 180s 两级提升至 600s） | 返回结果字典：包含 query_length/hits 键，hits 为列表（实际命中数取决于 NCBI 数据库，允许为空列表） |
| TC06-005 | 真实 XML 报文解析 | 轮询至 COMPLETE 取原始 XML 后 parse_result；未完成或 FAILED 自动重提（最多 3 次，间隔 15s，单次等待 120s——NCBI 偶发对连续提交限流，2026-10-15 登记） | 解析不抛异常；hits 中命中的 accession/identity_percent/e_value 字段存在 |
| TC06-006 | 真实别名解析 | resolve_official_symbol("ALK-7")（经 NCBI Gene 库） | 返回官方符号 ACVR1C（别名自动解析功能在线验证，2026-10-15 新增） |

## 已知兼容性发现（行为登记，不修改源码）

- NCBI 真实 XML 中查询信息标签为 `Iteration_query-def` / `Iteration_query-len`。`parse_result` 已于 2026-10-15 增加回退兼容（见 spec/changes/2026-10-15-P0问题修复.md），真实报文的 query_id / query_length 现可正常解析，Mock 层面验证见 TC04-011。

## 调试修复记录

- 在线用例首轮执行发现 `check_status` 缺陷：NCBI 任务完成时直接返回 XML 结果（无 `Status=COMPLETE` 标记），原实现永远无法返回 COMPLETE，导致 `run_blast` 真实访问必然超时。已按保留原逻辑方式修复（新增"响应为 XML 即视为完成"分支），详见 `spec/debug/2026-10-12-NCBI任务完成状态判断.md`。
