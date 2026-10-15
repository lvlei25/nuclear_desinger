# TC-05 综合评估与 GUI 集成用例

公共前置：S1 = `AUAGCUGCAUGAAAGUGCGCC`，S3 = `UUUUUUGGGGGCCCCCCCCCA`；GUI 与 demo 用例通过子进程执行，避免阻塞测试进程。

| 编号 | 用例标题 | 输入 | 预期结果 |
|------|---------|------|---------|
| TC05-001 | 综合评估-默认区域 | comprehensive_sirna_evaluation(S1) | final_score 约 76.37、recommendation='推荐'、target_region='起始密码子附近'、bonus=-5 |
| TC05-002 | 综合评估-指定 3'-UTR | target_region_name="3'-UTR" | final_score 约 91.37、bonus=10、recommendation='推荐' |
| TC05-003 | 综合评估-不推荐 | S3 | passes_hard_filters=False、recommendation='不推荐' |
| TC05-004 | 反义链一致性 | S1 评估结果 | antisense == reverse_complement(S1) == "GGCGCACUUUCAUGCAGCUAU" |
| TC05-005 | demo 命令行冒烟 | 子进程运行 `python demo.py` | 退出码 0，60 秒内完成 |
| TC05-006 | GUI 模块导入冒烟 | 子进程执行 `import gui` | 无异常退出，输出 GUI_IMPORT_OK（主窗口仅在 __main__ 启动） |
| TC05-007 | 端到端流程 | 设计 siRNA -> 综合评估 -> GalNAc 修饰方案 -> 结果存 CSV -> 读回 | 全链路无异常，CSV 读回包含 final_score 与修饰输出字段 |
| TC05-008 | 别名解析单元 | resolve_official_symbol("ALK-7")（Mock Gene 库 ESearch/ESummary） | 返回官方符号 "ACVR1C"（2026-10-15 新增） |
| TC05-009 | 检索回退流程 | search_nucleotide_mrna("ALK-7")（Mock 5 段响应：无命中→解析→重试→取序列） | 返回 1 条序列，重试请求 URL 使用官方符号 ACVR1C（2026-10-15 新增） |
| TC05-010 | GUI 源码静态审计 | 扫描 gui.py 源码 | 不存在未定义的 self.gene_name_var 引用，基因取值使用 gene_id_entry（2026-10-15 缺陷回归防护） |
