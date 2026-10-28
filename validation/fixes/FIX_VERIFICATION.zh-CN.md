# 科学一致性审计与修复复测

[English](FIX_VERIFICATION.md) · 审计日期：2026-10-28

**审计、算法与输入修复、回归测试及报告由 OpenAI Codex（ChatGPT）完成。**

结果：**PASS（仅针对下述已执行覆盖范围）**。这是用户授权修改算法后的新复测；[原审计](../audit/SCIENTIFIC_PARITY_AUDIT.md)和其结果未覆盖，原始 FAIL 结论作为修复前证据保留。

## 审计发现与判定

修复前，原 12 项 Rust 单元测试和原合成数据 9 项比较通过，但强化后的参数矩阵为
45 PASS / 7 FAIL，Student P 值网格有 124/252 项超出容差。原比较器缺文件、非有限值
和 ID 对齐检查不足，已增加必需输出、唯一标识、重排及负控制。

12K 原差异基因为 `ENSG00000149541`：R 为 turquoise，Rust 为 yellow。两边原始切树和
kME 清理标签一致。R 使用清理前的 eigengene 做重分配，Rust 原实现重新计算清理后
的 eigengene。相同基因的 Fisher 比值分别为 `1.0180061545e-6` 与 `9.1835531036e-7`，
位于 `1e-6` 阈值两侧；这不是已证实的固有浮点边界。修复遵循 R 的 eigengene 阶段，
没有更改阈值。[逐基因证据](../audit/artifacts/diagnostic_12000_comparison/mismatch_investigation.tsv)
和[阈值复算](../audit/artifacts/diagnostic_12000_comparison/reassignment_thresholds.tsv)保留供复核。

## 修复内容

1. `src/stats.rs`：不完全 beta 在收敛较好的区间使用继续分式，必要时使用互补变换和 log1p；弱相关区直接计算 r² 的补概率，避免在 1-r² 中丢失小量；达到迭代上限返回 NaN，而不是未收敛的概率。生产相关 P 值不再把极弱相关错误判为显著；正态 CDF 显式处理正负无穷，保证完美相关的 Fisher P 值为 0。
2. `src/network.rs`、`adjacency.rs`、`tom.rs`：整数幂仍用 powi，合法小数幂用 powf，不再截断；Nowick2 保留所选网络的 adjacency 变换，与实际 R WGCNA 1.74 行为一致。
3. `src/pipeline.rs`：清理和重分配共享原始动态模块的 eigengene；清理后不重新计算 eigengene，删除模块不作为候选。重分配后过小的来源模块按参考流程删除。原 0.3/0.5 kME 和 1e-6 重分配阈值未改。
4. `src/io.rs`、CLI、worker 和 benchmark：性状行按表达矩阵样本 ID 重排；重复、空白、缺失、意外 ID、矩阵维度或字段数不符显式报错。矩阵级计算 API 仍要求调用方提供已对齐的行。
5. `examples/real_bench.rs` 和测试副本：连续输出保留 17 位小数，避免六位小数文件无法验证 1e-9 的问题。没有提高比较容差。

新增两项重分配回归、三项样本对齐/坏输入回归；保留并修复原审计四项数值回归，另补充完美相关的 Student/Fisher P 值和正态 CDF 无穷端点回归。生产默认科学参数未修改，没有引入新的数值依赖。

## 实际执行结果

| group | PASS | FAIL | SKIPPED |
| --- | --- | --- | --- |
| cargo tests | 22 | 0 | 0 |
| parameter configurations | 52 | 0 | 0 |
| validator controls | 15 | 0 | 0 |
| all exported comparisons | 664 | 0 | 0 |
| real final module comparisons | 2 | 0 | 0 |


`cargo test --workspace` 的 14 项库测试、5 项数值/端点回归、3 项对齐测试全部通过。原先 45/52 的参数矩阵现在逐组重新生成 Rust/R 输出，不把旧 FAIL 改标签当成通过。原 seed-42 fixture 与确定性 challenging fixture 均沿用原审计，未更改生物数据。

| metric | status | total | max_abs | mean_abs | rmse | above_tolerance |
| --- | --- | --- | --- | --- | --- | --- |
| student | PASS | 252 | 2.30594709993426e-12 | 2.18335486266528e-13 | 4.80014529416613e-13 | 0 |
| fisher | PASS | 252 | 1.73107361778335e-12 | 1.8236716256901e-13 | 4.33714943321098e-13 | 0 |


Student P 值 252 个测试点由原先 124 FAIL 变为全部通过；Fisher 252 点也全部通过。绝对容差始终 1e-9，相对误差只作诊断。倒序性状表的新 Rust 输出与按 ID 对齐的 R 参考全部 9 项通过：[结果](artifacts/trait_reordered/comparison/parity_details.csv)。15 个验证器控制保持通过，故缺文件/错位/异常值仍会触发失败。

参数逐项结果：[parameter_coverage.csv](parameter_coverage.csv)。全部连续比较与错误统计：[parity_details.csv](parity_details.csv)。精确命令、耗时、退出码：[test_manifest.json](test_manifest.json)。

## 8K / 12K 真实数据

| case | total | exact_agreement | disagreeing | ARI | gene_order_identical | status |
| --- | --- | --- | --- | --- | --- | --- |
| diagnostic_12000_comparison | 12000 | 1 | 0 | 1 | TRUE | PASS |
| diagnostic_8000_comparison | 8000 | 1 | 0 | 1 | TRUE | PASS |


原 Demo FPKM/trait 文件、log10(x+1)、MAD top-N、power 6、unsigned/signed TOM/min、deepSplit 2、minModuleSize 30、merge .25、detectCutHeight .995 全部沿用。R 在独立进程中使用此前已确认可加载的 Accelerate BLAS（1 线程上限），Rust 使用默认可用线程。本轮是准确性诊断，不用这些耗时更新或宣称性能加速比。

| case | metric | status | max_abs | above_tolerance | total |
| --- | --- | --- | --- | --- | --- |
| artifacts/diagnostic_12000_comparison | module_trait_cor.tsv | PASS | 3.7470027081099e-15 | 0 | 55 |
| artifacts/diagnostic_12000_comparison | module_trait_p.tsv | PASS | 2.50000020685093e-12 | 0 | 55 |
| artifacts/diagnostic_12000_comparison | cor_probe.tsv | PASS | 9.753309271332e-14 | 0 | 780000 |
| artifacts/diagnostic_12000_comparison | tom_probe.tsv | PASS | 4.26880752968373e-14 | 0 | 780000 |
| artifacts/diagnostic_12000_comparison | sft.tsv | PASS | 7.32143234927207e-11 | 0 | 90 |
| artifacts/diagnostic_12000_comparison | tree.tsv | PASS | 6.56141807553468e-14 | 0 | 35997 |
| artifacts/diagnostic_12000_comparison | dynamic.tsv | PASS | 0 | 0 | 12000 |
| artifacts/diagnostic_12000_comparison | eigengenes.tsv | PASS | 3.60822483003176e-15 | 0 | 616 |
| artifacts/diagnostic_12000_comparison | kme.tsv | PASS | 1.00919272938427e-13 | 0 | 132000 |
| artifacts/diagnostic_12000_comparison | intramodular.tsv | PASS | 9.12336872715969e-11 | 0 | 12000 |
| artifacts/diagnostic_8000_comparison | module_trait_cor.tsv | PASS | 8.3544282603043e-15 | 0 | 30 |
| artifacts/diagnostic_8000_comparison | module_trait_p.tsv | PASS | 2.45051201552826e-12 | 0 | 30 |
| artifacts/diagnostic_8000_comparison | cor_probe.tsv | PASS | 5.42899059041702e-14 | 0 | 512000 |
| artifacts/diagnostic_8000_comparison | tom_probe.tsv | PASS | 4.71567229709535e-14 | 0 | 512000 |
| artifacts/diagnostic_8000_comparison | sft.tsv | PASS | 5.04769559483975e-11 | 0 | 90 |
| artifacts/diagnostic_8000_comparison | tree.tsv | PASS | 4.18554080283684e-14 | 0 | 23997 |
| artifacts/diagnostic_8000_comparison | dynamic.tsv | PASS | 0 | 0 | 8000 |
| artifacts/diagnostic_8000_comparison | eigengenes.tsv | PASS | 5.72979164115139e-15 | 0 | 336 |
| artifacts/diagnostic_8000_comparison | kme.tsv | PASS | 5.07371922253697e-14 | 0 | 48000 |
| artifacts/diagnostic_8000_comparison | intramodular.tsv | PASS | 6.10498318565078e-11 | 0 | 8000 |


12K 原差异基因 ENSG00000149541 的最终标签和参考 R 一并核对，详细文件见 `artifacts/diagnostic_12000_*/module_colors.txt` 和比较目录的不一致基因列表。本次通过不依赖放宽阈值、任意翻转 eigengene 或只按模块名称忽略差异。

## 性能审计方法

修复前审计区分默认 reference BLAS（A）与受控 Accelerate BLAS（B）。A 中 R WGCNA
为默认 1 线程，Rust 使用默认可用并行度；部分重复与审计任务重叠，只作为描述性结果。
B 中 R WGCNA/BLAS 和 Rust faer/Rayon 线程池各配置为 8；这不是对 OS 总线程数的断言。
非重叠计时和仅包含读取/归一化、MAD、完整模块检测和模块-性状关联，不重复加入单独
测量的 cor/TOM，也不包含输出写盘。默认 8K 三轮与默认 12K 一轮完成；12K 单轮约
834 秒，额外两轮显式标为 SKIPPED。B 在两个规模各重复三轮，记录中位数、范围与标准差。
这些是修复前代码的环境相关结果，不是修复后性能验证，也不是纯算法加速比。

详见[汇总](../audit/benchmark_summary.csv)、[逐轮记录](../audit/benchmark_details.csv)与
修复前审计的内存/负载限制。默认环境日志报告进程峰值 RSS；B 的每秒采样只提供下界。

## 覆盖限制与兼容性

- 大数据 cor/TOM 仍是约 64 行对全部基因的探针比较（8K 512000、12K 780000 元素），不宣称完整平方矩阵逐元素验证；全部基因的模块颜色、kME、模块内 connectivity 和树输出均比较。
- 真实数据 R 的 grey eigengene 在比较时显式排除，因为原 Rust 基准只导出非 grey；其他行列和基因身份严格检查。
- 本轮未重跑耗时较长的 default reference-BLAS 12K，也未重做性能重复实验；原审计中的此类证据仍单独保留。
- 此前无效的性状输入（重复 ID、缺行、额外列、依赖零填充等）现在会报错，是有意的输入契约收紧。表达/性状身份不匹配不能静默计算。
- 编译、测试与 git diff --check 通过。运行环境、工具与原有警告见命令日志。
- 这些结果支持列明范围内的修复，不代表所有参数/数据/缺失模式已证明等价。

## 来源与交付

基础提交：4f5055425ec9999dfe7b15422d24b5bb72cb1f7e；R 4.5.1 / WGCNA 1.74；Rust 1.99.1；Apple M2/16 GiB，与原审计同一机器。被测实现由 [source_hashes.json](source_hashes.json) 标识，文档与发布整理不改变被测数值算法。

原 `Demo/`、`bench_out/`、Cargo 依赖/锁文件及修复前 `validation/audit/` 证据均保留。[source_hashes.json](source_hashes.json) 标识生产源文件、测试和输入；[production_diff.patch](production_diff.patch) 列出代码变化。README 与 Quarto 区分修复前历史审计和修复后新结果；新的 HTML 单独输出为 `report/benchmark-fixed.html`。

复现命令与已发布/本地证据范围见 [validation/README.md](../README.md)。原始大矩阵、重复导出和本机符号链接保留在本地；Git 中提供汇总、命令日志、fixture、关键阶段数据与校验值。没有把未发布的大矩阵当成完整矩阵逐元素验证。
