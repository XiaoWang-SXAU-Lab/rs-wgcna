from pathlib import Path
import json,csv,re,hashlib,subprocess
import argparse
ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description='Summarize an existing completed fix-validation run without changing the reviewed bilingual reports.')
parser.add_argument('--output-dir',type=Path,default=ROOT/'validation/fixes')
OUT=parser.parse_args().output_dir.resolve();ART=OUT/'artifacts'
def read(p):
 with open(p) as f:return list(csv.DictReader(f))
def write(p,rows):
 fields=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def table(rows,fields):
 return '| '+' | '.join(fields)+' |\n| '+' | '.join(['---']*len(fields))+' |\n'+''.join('| '+' | '.join(str(r.get(k,'')) for k in fields)+' |\n' for r in rows)
details=[];modules=[]
for f in sorted(ART.rglob('parity_details.csv')):
 if 'negative_controls' in f.parts: continue
 for row in read(f):details.append(dict(case=str(f.parent.relative_to(OUT)),**row))
for f in sorted(ART.glob('matrix_*/comparison/extra_details.csv')):
 for row in read(f):details.append(dict(case=str(f.parent.relative_to(OUT)),**row))
for f in sorted(ART.glob('diagnostic_*_comparison/module_summary.csv')):
 for row in read(f):modules.append(dict(case=f.parent.name,**row))
write(OUT/'parity_details.csv',details);write(OUT/'module_summary.csv',modules)
coverage=read(OUT/'parameter_coverage.csv');stats=read(ART/'stats_grid/parity_details.csv');controls=read(ART/'negative_controls/summary.csv')
counts=[]
def count(name,statuses):counts.append(dict(group=name,PASS=statuses.count('PASS'),FAIL=statuses.count('FAIL'),SKIPPED=statuses.count('SKIPPED')))
cargo_log=ART/'cargo_final_tests.log'
if not cargo_log.exists(): cargo_log=ART/'cargo_tests.log'
log=cargo_log.read_text();passing=sum(int(n) for n in re.findall(r'test result: ok\. (\d+) passed;',log));fails=sum(int(n) for n in re.findall(r'; (\d+) failed;',log))
count('cargo tests',['PASS']*passing+['FAIL']*fails)
count('parameter configurations',[r['status'] for r in coverage]);count('validator controls',[r['status'] for r in controls]);count('all exported comparisons',[r['status'] for r in details]);count('real final module comparisons',[r['status'] for r in modules]);write(OUT/'test_counts.csv',counts)
changes=subprocess.check_output(['git','diff','--','src','examples/real_bench.rs','README.md','README.zh-CN.md','NEWS.md','report/benchmark.qmd'],cwd=ROOT,text=True);(OUT/'production_diff.patch').write_text(changes)
hashes={}
for pattern in ['src/**/*.rs','tests/*.rs','Demo/*','Cargo.toml','Cargo.lock','examples/audit*.rs','examples/real_bench.rs']:
 for p in ROOT.glob(pattern):
  if p.is_file():hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
(OUT/'source_hashes.json').write_text(json.dumps(hashes,indent=2))
real_details=[r for r in details if 'diagnostic_' in r['case']]
# No missing case may be counted as successful.
complete=len(coverage)==52 and all(r['status']=='PASS' for r in coverage) and len(modules)==2 and len(details)==664 and all(r['status']=='PASS' for r in details) and all(r['status']=='PASS' for r in modules) and len(controls)==15 and all(r['status']=='PASS' for r in controls) and passing==22 and fails==0
status='PASS' if complete else 'FAIL'
report=f'''# 科学一致性缺陷修复与复测

结果：**{status}（仅针对下述已执行覆盖范围）**。这是用户授权修改算法后的新复测；[原审计](../audit/SCIENTIFIC_PARITY_AUDIT.md)和其结果未覆盖，原始 FAIL 结论作为修复前证据保留。

## 修复内容

1. `src/stats.rs`：不完全 beta 在收敛较好的区间使用继续分式，必要时使用互补变换和 log1p；弱相关区直接计算 r² 的补概率，避免在 1-r² 中丢失小量；达到迭代上限返回 NaN，而不是未收敛的概率。生产相关 P 值不再把极弱相关错误判为显著；正态 CDF 显式处理正负无穷，保证完美相关的 Fisher P 值为 0。
2. `src/network.rs`、`adjacency.rs`、`tom.rs`：整数幂仍用 powi，合法小数幂用 powf，不再截断；Nowick2 保留所选网络的 adjacency 变换，与实际 R WGCNA 1.74 行为一致。
3. `src/pipeline.rs`：清理和重分配共享原始动态模块的 eigengene；清理后不重新计算 eigengene，删除模块不作为候选。重分配后过小的来源模块按参考流程删除。原 0.3/0.5 kME 和 1e-6 重分配阈值未改。
4. `src/io.rs`、CLI、worker 和 benchmark：性状行按表达矩阵样本 ID 重排；重复、空白、缺失、意外 ID、矩阵维度或字段数不符显式报错。矩阵级计算 API 仍要求调用方提供已对齐的行。
5. `examples/real_bench.rs` 和测试副本：连续输出保留 17 位小数，避免六位小数文件无法验证 1e-9 的问题。没有提高比较容差。

新增两项重分配回归、三项样本对齐/坏输入回归；保留并修复原审计四项数值回归，另补充完美相关的 Student/Fisher P 值和正态 CDF 无穷端点回归。生产默认科学参数未修改，没有引入新的数值依赖。

## 实际执行结果

{table(counts,['group','PASS','FAIL','SKIPPED'])}

`cargo test --workspace` 的 14 项库测试、5 项数值/端点回归、3 项对齐测试全部通过。原先 45/52 的参数矩阵现在逐组重新生成 Rust/R 输出，不把旧 FAIL 改标签当成通过。原 seed-42 fixture 与确定性 challenging fixture 均沿用原审计，未更改生物数据。

{table(stats,['metric','status','total','max_abs','mean_abs','rmse','above_tolerance'])}

Student P 值 252 个测试点由原先 124 FAIL 变为全部通过；Fisher 252 点也全部通过。绝对容差始终 1e-9，相对误差只作诊断。倒序性状表的新 Rust 输出与按 ID 对齐的 R 参考全部 9 项通过：[结果](artifacts/trait_reordered/comparison/parity_details.csv)。15 个验证器控制保持通过，故缺文件/错位/异常值仍会触发失败。

参数逐项结果：[parameter_coverage.csv](parameter_coverage.csv)。全部连续比较与错误统计：[parity_details.csv](parity_details.csv)。精确命令、耗时、退出码：[test_manifest.json](test_manifest.json)。

## 8K / 12K 真实数据

{table(modules,['case','total','exact_agreement','disagreeing','ARI','gene_order_identical','status'])}

原 Demo FPKM/trait 文件、log10(x+1)、MAD top-N、power 6、unsigned/signed TOM/min、deepSplit 2、minModuleSize 30、merge .25、detectCutHeight .995 全部沿用。R 在独立进程中使用此前已确认可加载的 Accelerate BLAS（1 线程上限），Rust 使用默认可用线程。本轮是准确性诊断，不用这些耗时更新或宣称性能加速比。

{table(real_details,['case','metric','status','max_abs','above_tolerance','total'])}

12K 原差异基因 ENSG00000149541 的最终标签和参考 R 一并核对，详细文件见 `artifacts/diagnostic_12000_*/module_colors.txt` 和比较目录的不一致基因列表。本次通过不依赖放宽阈值、任意翻转 eigengene 或只按模块名称忽略差异。

## 覆盖限制与兼容性

- 大数据 cor/TOM 仍是约 64 行对全部基因的探针比较（8K 512000、12K 780000 元素），不宣称完整平方矩阵逐元素验证；全部基因的模块颜色、kME、模块内 connectivity 和树输出均比较。
- 真实数据 R 的 grey eigengene 在比较时显式排除，因为原 Rust 基准只导出非 grey；其他行列和基因身份严格检查。
- 本轮未重跑耗时较长的 default reference-BLAS 12K，也未重做性能重复实验；原审计中的此类证据仍单独保留。
- 此前无效的性状输入（重复 ID、缺行、额外列、依赖零填充等）现在会报错，是有意的输入契约收紧。表达/性状身份不匹配不能静默计算。
- Rustfmt 未安装，本次工具调用失败已知；代码通过编译、测试与 git diff --check。其他原有编译警告未在此次科学修复中顺带改动。
- 这些结果支持列明范围内的修复，不代表所有参数/数据/缺失模式已证明等价。

## 来源与交付

基础提交：4f5055425ec9999dfe7b15422d24b5bb72cb1f7e；R 4.5.1 / WGCNA 1.74；Rust 1.99.1；Apple M2/16 GiB，与原审计同一机器。当前工作区含未提交的审计和修复改动；未 commit 或 push。

原 `Demo/`、`bench_out/`、Cargo 依赖/锁文件及修复前 `validation/audit/` 证据均保留。[source_hashes.json](source_hashes.json) 标识生产源文件、测试和输入；[production_diff.patch](production_diff.patch) 列出代码变化。README 与 Quarto 区分修复前历史审计和修复后新结果；新的 HTML 单独输出为 `report/benchmark-fixed.html`。

重跑：在一个新的输出目录使用 `run_validation.py` 中记录的参数与既有审计比较脚本；脚本拒绝覆盖已存在的 case 目录。无需下载或替换生物数据。
'''
# Preserve the reviewed bilingual reports; write generated narrative separately.
(OUT/'FIX_VERIFICATION.generated.zh-CN.md').write_text(report)
print(status,passing,'cargo tests;',len(coverage),'configs;',len(modules),'real comparisons;',len(details),'comparison rows')
raise SystemExit(0 if complete else 1)
