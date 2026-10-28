from run_audit import ROOT,A
from aggregate_results import read
import pathlib,subprocess,json,csv

def table(rows,cols):
 if not rows:return '（暂无已完成结果；不得视为 PASS。）\n'
 def fmt(v):
  if isinstance(v,float):return f'{v:.10g}'
  return str(v).replace('|','/').replace('\n',' ')
 return '| '+' | '.join(cols)+' |\n| '+' | '.join(['---']*len(cols))+' |\n'+''.join('| '+' | '.join(fmt(r.get(k,'')) for k in cols)+' |\n' for r in rows)
counts=read(A/'test_counts.csv');coverage=read(A/'parameter_coverage.csv');summary=read(A/'parity_summary.csv');details=read(A/'parity_details.csv');bench=read(A/'benchmark_summary.csv')
real=[r for r in summary if r['kind']=='module' and any(k in r['case'] for k in ['real_','perfB_','diagnostic_'])]
core=[]
for metric in ['cor.tsv','adjacency.tsv','tom.tsv','sft.tsv','eigengenes.tsv','kme.tsv','module_trait_cor.tsv','module_trait_p.tsv']:
 vals=[float(r['max_abs']) for r in details if r['metric']==metric and ('matrix_' in r['case'] or 'hardened_final' in r['case']) and r['max_abs'] not in ('NA','')]
 core.append(dict(metric=metric,max_abs=max(vals),failed_comparisons=sum(r['status']=='FAIL' for r in details if r['metric']==metric and ('matrix_' in r['case'] or 'hardened_final' in r['case']))))
large=[r for r in details if 'diagnostic_' in r['case']]
stats=read(A/'artifacts/stats_grid/parity_details.csv')
stages=read(A/'stage_summary.csv')
final_log=A/'artifacts/cargo_verified_final.log'; final_regressions=4 if '4 failed' in final_log.read_text() else '见实际日志'
q='PASS' if (A/'artifacts/quarto_render_successful.execution.json').exists() and json.loads((A/'artifacts/quarto_render_successful.execution.json').read_text())['exit_code']==0 else 'NOT EXECUTED'
text=f'''# rs-WGCNA 科学一致性审计

审计日期：2026-10-28（Asia/Shanghai）。最终科学判定：**FAIL**。审计过程保留全部失败，不改变科学容差或 Rust 生产算法。

## 1. 执行摘要

本地 `rs_WGCNA` 在任务开始时与 GitHub `master`/远端 HEAD 完全一致，提交为 `4f5055425ec9999dfe7b15422d24b5bb72cb1f7e`，初始工作区干净。依据：[environment.log](artifacts/environment.log)。这次新增的是验证脚本、测试和报告，未提交或推送。结束时的远端复查遇到 HTTP2 网络错误，HTTP1.1 重试在 20 秒超时，标为 SKIPPED；一致性结论有任务开始时的成功远端证据，不冒充结束时重新验证通过。见 [remote_final.json](remote_final.json)。

原始 12 个 Rust 单元测试及原始合成数据 9 项比较通过；但原比较器存在会掩盖缺文件、非有限值、ID 错位的漏洞。强化后，52 组参数配置为 **{sum(r['status']=='PASS' for r in coverage)} PASS / {sum(r['status']=='FAIL' for r in coverage)} FAIL**。默认整数幂、常用 unsigned/signed TOM 的小型矩阵通常达到很高精度，但不能据此推断所有参数及统计函数等价。

主要发现：小数幂被截断；`signed Nowick 2` 的 TOM 映射不一致；Student P 值在弱相关区有严重数值误差；带 ID 的性状表重排会导致 Rust 入口按错误样本配对；12K 的单基因差异与清理后重新计算 eigengene 的重分配路径有关。8K 最终模块为 100% 一致，12K 为 11999/12000（99.9916667%），ARI 约 0.999808868。现有证据不支持 README 的“仅为固有浮点边界、不是逻辑问题”断言。

核心证据：[parity_details.csv](parity_details.csv)、[parity_summary.csv](parity_summary.csv)、[test_counts.csv](test_counts.csv)、[test_manifest.json](test_manifest.json)。原始基准文件未覆盖。

## 2. 精确源代码与运行环境

- GitHub：`https://github.com/XiaoWang-SXAU-Lab/rs-wgcna`；工作目录：`{ROOT}`。
- Rust/Cargo 1.99.1，release 编译使用原 Cargo 配置和锁文件；测试使用 cargo test 的默认 test profile。
- R 4.5.1；WGCNA 1.74；dynamicTreeCut 1.63.1；fastcluster 1.3.0。WGCNA 安装包提示 built under R 4.5.2，这一版本差异已保留在日志。
- Apple M2，8 逻辑核，16 GiB RAM，macOS Tahoe 26.5，aarch64。
- Scenario A：R 默认 libRblas（reference BLAS），WGCNA 1 线程；Rust 按可用并行度配置 faer，8 线程，原入口没有显式约束全局 Rayon 池。
- Scenario B：进程级加载已安装 Accelerate BLAS，WGCNA 与 BLAS 各配置 8 线程；测试专用 Rust 基准把 faer 和 Rayon 均配置为 8。这里是各线程池配置，不是对总 OS 线程数的断言。
- 真实数据诊断 R 使用 Accelerate，BLAS 上限 1、WGCNA 1。它与默认 reference BLAS 重跑的最终模块结果分开报告。

完整信息：[environment.json](environment.json)、[R_environment_default.json](artifacts/R_environment_default.json)、[R_environment_accelerate8.json](artifacts/R_environment_accelerate8.json)。[file_hashes.json](file_hashes.json) 保存生产源码、输入与历史文件的 SHA-256。未下载替代生物数据。

## 3. 原测试基线与实际执行

修改验证脚本之前已运行 `cargo test --workspace`：12 PASS、0 FAIL、0 ignored。原脚本 `run_parity.sh 6 1e-9` 在隔离输出目录中运行，9 项比较 PASS，最大误差 1.59394719645434e-12（模块–性状 P 值）。原脚本副本和输出保存在 [artifacts/baseline](artifacts/baseline)。初始基础设施审计写于更改脚本之前：[INITIAL_AUDIT.md](INITIAL_AUDIT.md)。

强化后重新运行相同核心流程：[parity_hardened_final.log](artifacts/parity_hardened_final.log)，退出码 0。最终 `cargo test --workspace` 中原有 12 项仍通过，新增 {final_regressions} 项未忽略的回归测试按预期暴露缺陷，命令非零退出；这是真实 FAIL，不记为 PASS。

{table(counts,['group','PASS','FAIL','SKIPPED','NOT SUPPORTED','NOT EXECUTED'])}

计数按表中层级分别统计：参数配置、比较项、网格元素、控制测试不能相加当成独立实验总数；历史输出不计入当前验证通过数。负控制的 PASS 表示错误被检测到并返回预期非零退出码。

精确命令、目录、时长、退出码：[test_manifest.json](test_manifest.json)；完整输出：[test_execution_log.txt](test_execution_log.txt)。过程中发现并修正的测试辅助程序解析/编译错误保留在相应日志，最终证据使用 verified/final 运行；不把失败尝试抹去。

## 4. 原比较器漏洞

`scripts/compare.R` 原来把必需文件缺失当成 SKIP，缺失颜色文件甚至无提示；只检查矩阵尺寸，不检查行列 ID；按位置相减；max/sum 使用 na.rm=TRUE；仅报 max absolute error；不报 ARI、分区列联表或重排策略。混合有限值/NA 时可隐藏错误。原模块标签向量没有基因 ID。原 Quarto 的动态切树 100% 与部分性能/科学结论是手工写入。

真实数据基准另有两个关键约束：Rust 只输出非 grey 模块，R 包括 grey；Rust 原基准的连续结果只输出六位小数，无法验证 1e-9 精度。两点均显式记录，未将不相同的文件直接按位置比较。

## 5. 验证改动与失败语义

强化比较器要求全部八个连续文件和颜色文件存在；验证唯一非空 ID、尺寸与标识集合，并按 ID 显式对齐。表达矩阵提供额外的基因/样本集合约束。无 ID 的颜色向量使用各自已标注的 correlation 行顺序，并核对原表达矩阵顺序。

所有非有限值（NA、NaN、Inf、-Inf）均失败；当前确定性 fixture 没有合理的预期缺失，因此没有放宽掩码策略。每个连续比较保存总数、max/mean absolute error、RMSE、最大相对误差、超容差数量/比例。绝对容差保持 1e-9；相对分母为 max(abs(R),1e-12)，近零相对误差不单独决定 FAIL。SFT 仅规范化已知列名别名，没有对未知列名宽松匹配。

颜色精确一致率与 ARI 分开；输出基因不一致列表和列联表。必需失败返回非零，缺证据不会变为 PASS。`run_parity.sh` 现在默认新建输出目录，保留旧 r_out/rust_out/testdata。测试扩展继续调用现有 Rust/R 核心，没有重新实现网络算法。

15 项控制均符合预期，包括缺文件、错误数值、全部四种非有限值、重复/错误 ID、错误样本 ID、标签错位、缺颜色，以及正确重排、手算误差与颜色任意重命名的正控制。[控制结果](artifacts/negative_controls/summary.csv)。

## 6. 参数覆盖

{table(coverage,['case','network','TOM','denominator','power','deepSplit','mergeCutHeight','fixture','status'])}

网络类型为 unsigned / signed / signed hybrid；六个非 none TOM 类型 × 两种分母均测试，另含 none。幂为 1、6、12 及实际被接口接受的小数 2.5。deepSplit 0/1/2/3/4，mergeCutHeight 0/.15/.25/.5；不是所有参数的笛卡尔积。原 seed-42 的 150×40 fixture、确定性弱信号/反相关扰动 fixture 用于参数覆盖；大模块由原 8K/12K 数据覆盖。

所有表内组合均被两边接受，失败不能标为 NOT SUPPORTED。相关类型仅测 Pearson；不宣称 bicor、分块、多缺失值、外部 branch split、cosine 等未覆盖能力已验证。`signed Nowick 2` 六组 FAIL；fractional power 一组 FAIL。部分失败的颜色名称不一致但 ARI=1，这是分区重命名而非生物分区改变。

## 7. 连续量误差与独立统计网格

以下为小型参数测试与强化核心流程的最大误差，**包括已知失败参数**，不能解读为通过配置的精度：

{table(core,['metric','max_abs','failed_comparisons'])}

独立统计网格：9 个样本量（3 到 1000），14 个相关幅度，正负两种符号，共 252 项（零值正负重复保留，不把它们称为独立随机样本）。

{table(stats,['metric','status','total','max_abs','mean_abs','rmse','above_tolerance'])}

Student grid 中 124/252 超过 1e-9；r=1e-8、n=5 时 R≈0.9999999873，Rust≈4.217e-6。Fisher 的绝对误差通过。源代码中的 incomplete-beta 继续分式缺少近 1 区域的互补变换/明确不收敛处理，是可疑数值路径；本审计只证明失败，不擅自修复。[PVALUE_MRE.md](PVALUE_MRE.md)。

小数幂最大 adjacency 误差约 0.07048239、TOM 0.08142745：[FRACTIONAL_POWER_MRE.md](FRACTIONAL_POWER_MRE.md)。Nowick2 TOM 最大误差约 0.23192776，2 基因例子产生 +0.5 / -0.5：[NOWICK2_MRE.md](NOWICK2_MRE.md)。性状重排使模块–性状相关最大误差约 0.69805869：[TRAIT_ALIGNMENT_MRE.md](TRAIT_ALIGNMENT_MRE.md)。

## 8. 精确标签与置换不变分区

{table(real,['case','total','exact_agreement','disagreeing','ARI','gene_order_identical','status'])}

所有真实数据比较先验证唯一基因 ID 和一致基因集，并按 ID 对齐。MAD 保留顺序也核对一致。输入包含 36953 基因、56 样本；原性状顺序与表达矩阵样本顺序完全一致：[input_alignment.log](artifacts/input_alignment.log)。因此该固定输入不受性状重排缺陷影响。基因/样本 ID 不匹配不能依靠矩阵数值看起来接近而通过。

## 9. 8K / 12K 端到端验证

原 Demo 数据、log10(x+1)、MAD top-N、power=6、unsigned network、signed TOM、min denominator、minModuleSize=30、merge=.25、deepSplit=2、detectCutHeight=.995，未为提高通过率调整。

默认环境 8K 完整运行 3 次，12K 完整运行 1 次；受控优化环境重复次数见 benchmark_summary。12K reference-BLAS 单轮实际约 834 秒，因此另外两轮默认 12K 标 SKIPPED，改用优化环境重复及诊断运行保留覆盖。[repetition_policy.json](repetition_policy.json)。初始运行驱动在第一轮结束后通过预留的输出保护目录停止，guard 的非零终止是执行控制，不是额外一次生物学 FAIL；未把预留目录当成结果。

诊断额外比较 SFT、cor/TOM 探针、全部树 merge/height、动态切树、清理、最终 eigengene/kME、模块–性状 r/p、模块内 connectivity。8K cor/TOM 各 512000 元素，12K 各 780000 元素；这是均匀选取约 64 行对全部基因的探针，**不是完整平方矩阵逐元素验证**。大矩阵 adjacency 没有导出，NOT EXECUTED；小型参数矩阵的 adjacency 完整比较。

{table(large,['case','metric','status','max_abs','above_tolerance','total'])}

8K 模块一致，但高精度 P 值仍有约 4.10e-8 误差；12K 模块不同会引发 eigengene/kME/connectivity 的后续差异，不能全部解释为线性代数舍入。原 benchmark 的六位小数输出达到不了 1e-9，严格比较仍如实 FAIL；没有通过提高容差把它们变绿。

## 10. 12K 差异基因与最早分歧

唯一重现的基因是 **ENSG00000149541**：R=turquoise，Rust=yellow。历史输出也有同一基因。原始动态标签与清理标签均为 1；R premerge 为 turquoise。

{table(stages,['genes','stage','total','exact_agreement','ARI','disagreeing','status'])}

R 使用原动态模块的 eigengene 执行 Fisher 比值重分配；Rust 的私有重分配函数重新计算清理后 eigengene。实际公共 API 复测与 R 独立重算一致：原动态模块 ratio≈1.0180061545e-6（保留），清理后 ratio≈9.1835531036e-7（重分配），生产阈值均为 1e-6。最佳 kME 从约 .72014109 变成 .72114484，差异约 .001，远大于浮点误差。两边同一 eigengene 基础上的 kME 仅有机器精度量级差别。

因此证据指向**清理与重分配之间的 eigengene 基础不同**，最早已证实的决策分歧为重分配；没有证据支持树 ties、标签任意重命名或 eigengene 方向是此基因分歧的主要原因。Rust 私有函数未直接导出 premerge 状态，故保留这一跟踪限制，建议后续增加专门可观测接口。

详情：[REASSIGNMENT_FINDING.md](REASSIGNMENT_FINDING.md)、[逐基因 kME](artifacts/diagnostic_12000_comparison/mismatch_investigation.tsv)、[R 阈值重算](artifacts/diagnostic_12000_comparison/reassignment_thresholds.tsv)、[Rust 实测](artifacts/diagnostic_12000_comparison/reassignment_rust_measurements_verified.tsv)。最终 own kME 距离 minKMEtoStay 的值已导出，但这个清理阈值不是引发当前分歧的 Fisher ratio 阈值，二者没有混用。

## 11. 性能方法与结果

计时边界沿用原 real_bench：读取与 log10、MAD 选基因、cor、TOM、完整模块检测、模块–性状。完整模块检测内部重算 cor/TOM；端到端和仅求 log10 + MAD + module_detection_total + module_trait，不重复计入单独 cor/TOM。它并非程序全部壁钟（输出写盘未计入）。两种实现的相同性状模块集不完全一致（grey），读取实现与线程池也不同，不能声称算法单独加速。

{table([r for r in bench if r['step'] in ['nonoverlap_total','tom','module_detection_total']],['scenario','genes','implementation','step','repeats','median_seconds','min_seconds','max_seconds','sd_seconds','threads','BLAS'])}

默认环境的部分重复与编译、参数测试或诊断并发，存在负载污染，均为描述性数据，不能用于可靠速度倍数。受控 B 的基准驱动等待默认测试与大规模诊断结束，逐轮串行；桌面其他活动不受控制，不是专用性能实验。优化 BLAS 已通过真实加载路径验证，未修改系统链接。

默认 `/usr/bin/time -l` 日志含进程壁钟及最大 RSS/peak footprint；不能将单进程 RSS 简单等同于所有子进程/整个系统内存。B 的 [perfB_sampled_memory.json](artifacts/perfB_sampled_memory.json) 为 1 秒间隔 RSS 采样下界，启动晚于首个 Rust repeat，不冒充精确峰值。原始历史 timings 保存在 bench_out，源提交/线程环境未知，仅作历史对照。

完整逐轮计时：[benchmark_details.csv](benchmark_details.csv)，中位数、范围和标准差：[benchmark_summary.csv](benchmark_summary.csv)。不从单次 12K 默认测量推断方差，不丢弃慢轮改善速度比。

## 12. FAIL / SKIPPED / NOT SUPPORTED / NOT EXECUTED

- FAIL：7 参数配置；Student grid；性状重排 MRE；12K 最终模块；真实连续量的明确误差或原输出精度不足；新增 4 个回归测试。具体每项和退出码见 CSV/manifest。
- SKIPPED：另外两次昂贵的默认 reference-BLAS 12K 重复，原因明确记录。
- NOT SUPPORTED：当前枚举测试中没有此类组合。没有把已接受但算错的参数改标 unsupported。
- NOT EXECUTED：大数据完整 adjacency/cor/TOM 平方矩阵逐元素比较；生产 Rust 私有函数的直接 premerge/reassignment 状态导出；原基准未导出的连续/中间量（在其原始运行上）。探针和单独诊断不能抹掉这些缺证据项。

## 13. 已知限制

这是有限的确定性配置与两个来自同一 Demo 源矩阵的基因子集验证，不代表所有生物数据、噪声/缺失模式或参数。8K 和 12K 不是独立来源数据。未测的可选功能不能当成通过。eigengene 方向未通过任意翻转“修正”：沿用各自默认平均表达方向并报告原始符号；ARI 不依赖颜色名。实际数据阶段诊断使用优化 BLAS，默认 BLAS 仅最终结果复核，因此中间阶段环境差异透明保留。

历史 README 中精度值准确重现于原合成 fixture，却不足以支持通用“机器精度”等价。旧 HTML 保留为历史证据，不能当作当前自动生成的报告。Quarto 新报告状态目前为 {q}，详见后述最终检查。

## 14. 最终科学判断与生产文件验证

**整体 FAIL；有限配置的核心相关/邻接/TOM 通过，不等于全实现科学等价。** Student P 值问题具有实际显著性判断风险；Nowick2、小数幂与重分配差异均有可执行证据。无需也没有通过调整容差、改数据或算法来隐藏失败。

`git diff --exit-code -- src Cargo.toml Cargo.lock Demo bench_out report/benchmark.html` 返回 0：[unchanged_production.log](artifacts/unchanged_production.log)。完整变更：[final_git_diff.patch](final_git_diff.patch)；状态：[final_git_status.txt](final_git_status.txt)。新增例子是测试专用 harness，回归测试未忽略；源算法未修改，未 commit/push。

`scripts/compare.R`、`scripts/r_reference.R`、`scripts/run_parity.sh` 和 `report/benchmark.qmd` 是原有文件改动；新增 tests/audit_numeric_contract.rs、examples/audit*.rs、scripts/audit_real_reference.R 和 validation/audit 下的运行/统计/证据文件。Quarto 读取真实 CSV，不再手写动态切树或整体等价结论；生成到 report/benchmark-audit.html，旧 report/benchmark.html 未覆盖。渲染状态：{q}。

## 15. 建议后续动作

1. 单独修复并严格回归测试 Student P 值，覆盖近零相关、端点和多样本量；对不收敛显式报告。
2. 明确 fractional power 支持契约，支持实数幂或明确拒绝；核对 Nowick2 的 WGCNA C 端映射与各种网络类型。
3. 使重分配时的 eigengene 生成阶段与参考实现一致，并加入 ENSG00000149541 回归 fixture；不要修改阈值使该基因“通过”。
4. Rust 输入入口按样本 ID 对齐性状，或明确拒绝顺序不符；禁止静默丢弃身份信息。
5. 完善公开阶段诊断、完整大矩阵流式比较及独立真实数据集验证，再更新 README 的科学/性能结论。当前审计未改 README，避免在验证任务中混入生产宣传改动。

复现说明：[REPRODUCE.md](REPRODUCE.md)。审计数据包目录即 `validation/audit/`；所有关键结论可沿相对链接追溯到机器可读数据或原始可执行命令。
'''
(A/'SCIENTIFIC_PARITY_AUDIT.md').write_text(text)
print('report generated',len(text),'characters')
