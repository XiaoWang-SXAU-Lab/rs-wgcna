# Scientific parity audit and fix verification

[中文](FIX_VERIFICATION.zh-CN.md) · Audit date: 2026-10-28 (Asia/Shanghai)

**Scientific auditing, implementation fixes, regression tests and reporting were performed by OpenAI Codex (ChatGPT), at the maintainer's request.**

## Assessment

**PASS within the executed coverage below.** The pre-fix audit remains a separate
[historical FAIL assessment](../audit/SCIENTIFIC_PARITY_AUDIT.md). No failed or missing
historical comparison has been relabelled as a passing post-fix result. The numerical
algorithms and inputs were corrected, and the affected comparisons were executed again.

The reference is R WGCNA 1.74. Both real-data subsets have exact final colour-label
agreement and ARI=1, including all 12,000 genes in the larger subset. This is finite,
fixture-specific evidence, not a proof of equivalence for every dataset or parameter.

## Source and environment

- Baseline source: `4f5055425ec9999dfe7b15422d24b5bb72cb1f7e`.
- Tested implementation content: [source_hashes.json](source_hashes.json).
- Rust/Cargo 1.99.1; R 4.5.1; WGCNA 1.74; dynamicTreeCut 1.63.1;
  fastcluster 1.3.0. The installed WGCNA package reports that it was built under
  R 4.5.2; this warning is preserved in the logs.
- Apple M2, 8 logical cores, 16 GiB RAM, macOS Tahoe 26.5, arm64.
- Real-data accuracy diagnostics: R uses process-local Accelerate BLAS with a
  one-thread cap; Rust uses its default available parallelism. System R libraries
  were not changed. These diagnostic durations are not performance measurements.
- Cargo manifest/lockfile, original `Demo/` data and historical `bench_out/` remain
  unchanged. [Preservation check](artifacts/preserved_inputs.json) verifies 28 hashes.

## Audit findings and corrections

The original 12 Rust tests and nine synthetic comparisons passed. Strengthened
validation exposed seven failing configurations out of 52 and 124 failing Student
P-value cases out of 252. The original comparator could silently skip absent outputs,
ignore non-finite values or compare differently ordered identifiers positionally.
Required-file, identifier, non-finite and failure-exit checks now have explicit controls.

### Student and Fisher P-values

The former incomplete-beta evaluation could return an unconverged fraction for weak
correlations. For `r=1e-8, n=5`, R gives approximately `0.9999999873`, whereas the old
Rust result was approximately `0.0000042171`. The audited grid's maximum absolute
Student P-value error was `0.9999957702`.

`src/stats.rs` now uses the complementary incomplete-beta expansion in its favourable
convergence region, `log1p`, and a direct `r²` complement for weak correlations. A
non-converged fraction returns NaN rather than a spurious probability. Infinite Fisher
z statistics are handled explicitly, giving P=0 for perfect correlation. No scientific
significance, assignment or comparison threshold was relaxed.

### Fractional powers and Nowick2 TOM

Adjacency and TOM previously truncated an accepted real power to `i32`. At correlation
0.25 and power 2.5, the correct unsigned adjacency is 0.03125; truncation produced
0.0625. `Power` retains the integer `powi` fast path and uses `powf` for real powers.

For `signed Nowick 2`, Rust applied the original Nowick sign-preserving adjacency
override to the second form as well. In the two-gene anticorrelation reproducer,
R returns an off-diagonal TOM of +0.5 while the former Rust output was -0.5. The second
form now retains the chosen network's adjacency mapping, matching the observed R
WGCNA 1.74 reference across all six tested network/denominator configurations.

Executable reproductions are retained in `tests/audit_numeric_contract.rs`; the
pre-fix [fractional-power](../audit/FRACTIONAL_POWER_MRE.md),
[Nowick2](../audit/NOWICK2_MRE.md), and [P-value](../audit/PVALUE_MRE.md) records remain available.

### Sample identity

The CLI, worker and benchmark previously discarded uploaded trait sample IDs. Reversing
otherwise correctly labelled trait rows changed module-trait correlations by as much as
0.698 in the audited example. Trait rows are now aligned to expression sample IDs.
Duplicate/blank IDs, missing/unexpected samples, mismatched dimensions and malformed
row widths fail explicitly. Matrix-only numerical APIs continue to require aligned rows.

The reversed-trait fixture passes all nine comparisons:
[comparison results](artifacts/trait_reordered/comparison/parity_details.csv).

### The 12K reassignment discrepancy

The reproduced discordant gene was `ENSG00000149541`: R assigned turquoise, Rust yellow.
Raw dynamic-tree-cut and post-cleanup labels agreed. R reassignment reused the original
dynamic-module eigengenes; Rust recomputed eigengenes from cleaned modules.

| Eigengene basis | Own absolute kME | Best absolute kME | P(best)/P(own) | Decision at 1e-6 |
|---|---:|---:|---:|---|
| Original dynamic modules | 0.5124264102 | 0.7201410897 | 1.0180061545e-6 | retain |
| Cleaned modules | 0.5124214676 | 0.7211448429 | 9.1835531036e-7 | reassign |

The kME change is about 0.001, well above roundoff. Rust now reuses the original
module eigengenes, excludes deleted modules from candidates, and deletes source
modules made too small by reassignment. The 0.3/0.5 kME criteria and `1e-6`
reassignment threshold are unchanged. The corrected final label is turquoise;
all 12,000 final labels now match R exactly.

Evidence: [gene-level stages and kME](../audit/artifacts/diagnostic_12000_comparison/mismatch_investigation.tsv),
[R threshold recomputation](../audit/artifacts/diagnostic_12000_comparison/reassignment_thresholds.tsv),
[Rust measurements](../audit/artifacts/diagnostic_12000_comparison/reassignment_rust_measurements_verified.tsv),
and [post-fix disagreements](artifacts/diagnostic_12000_comparison/module_disagreements.tsv).

### Output precision

Benchmark continuous outputs now retain 17 decimal places. Six-decimal output cannot
establish `1e-9` agreement; the old quantized files remain historical evidence. The
comparison tolerance was not increased to accommodate them.

## Executed test results

| group | PASS | FAIL | SKIPPED |
| --- | --- | --- | --- |
| cargo tests | 22 | 0 | 0 |
| parameter configurations | 52 | 0 | 0 |
| validator controls | 15 | 0 | 0 |
| all exported comparisons | 664 | 0 | 0 |
| real final module comparisons | 2 | 0 | 0 |


Cargo results comprise 14 library tests, five numerical/endpoint regressions and three
trait-alignment/input tests. Fifteen validator controls include deliberately missing
files, perturbed values, NA/NaN/Inf/-Inf, duplicate/wrong IDs and shuffled labels, plus
positive checks for correct reordering, numerical error calculations and arbitrary
partition renaming. A passing negative control means the intended error was detected
and returned its expected nonzero exit code.

The 52 configurations were regenerated in Rust and R. Coverage includes unsigned,
signed and signed-hybrid networks; all six non-none TOM forms plus none; min/mean
denominators; powers 1, 6, 12 and 2.5; deepSplit 0–4; and representative merge heights.
It is a selected matrix, not their full Cartesian product. The deterministic 150×40
seed-42 fixture and the same challenging fixture were reused unchanged.

[Parameter coverage](parameter_coverage.csv) · [All error statistics](parity_details.csv) ·
[Commands and exit codes](test_manifest.json) · [Validator controls](artifacts/negative_controls/summary.csv).
Counts at different levels overlap and must not be summed as independent experiments.

## Numerical policy and statistical grid

Required files must exist. Dimensions, unique nonempty row/column IDs and identifier
sets are validated before explicit alignment. Unlabelled synthetic colour vectors
follow the input gene order, checked against labelled correlation output. Unexpected
NA, NaN or infinities fail. No `na.rm` exclusion hides numerical discrepancies.

The absolute tolerance remains `1e-9`. Each continuous comparison records maximum and
mean absolute error, RMSE, maximum relative error, count/proportion exceeding tolerance,
and total elements. Relative error uses `max(abs(reference),1e-12)` and is diagnostic;
near-zero relative error is not a separate failing criterion. Exact labels and
permutation-invariant ARI are reported separately, with contingency tables.

The statistical grid contains nine sample counts from 3 to 1000, 14 correlation
magnitudes and both signs (252 cases; positive/negative zero are retained duplicates).

| metric | status | total | max_abs | mean_abs | rmse | above_tolerance |
| --- | --- | --- | --- | --- | --- | --- |
| student | PASS | 252 | 2.30594709993426e-12 | 2.18335486266528e-13 | 4.80014529416613e-13 | 0 |
| fisher | PASS | 252 | 1.73107361778335e-12 | 1.8236716256901e-13 | 4.33714943321098e-13 | 0 |


Both 252-case grids pass; separate endpoint regressions verify r=0 and r=±1.

## Real-data protocol and results

Both implementations use the original 36,953-gene × 56-sample `Demo/` matrix and trait
table: `log10(x+1)`, MAD top 8K/12K, power 6, unsigned network, signed TOM, min
denominator, minimum module size 30, merge height 0.25, deepSplit 2 and detection
height 0.995. Gene/sample IDs are unique; the original trait table has identical
sample order. Selected gene identities and ordering are checked independently.

| case | total | exact_agreement | disagreeing | ARI | gene_order_identical | status |
| --- | --- | --- | --- | --- | --- | --- |
| diagnostic_12000_comparison | 12000 | 1 | 0 | 1 | TRUE | PASS |
| diagnostic_8000_comparison | 8000 | 1 | 0 | 1 | TRUE | PASS |


The exact-agreement value 1 means 100%. Every final gene label is checked, including
grey. R's additional grey eigengene is explicitly excluded only from numerical
module-level comparisons because the Rust benchmark exports non-grey eigengenes.

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


The complete tree merge/height exports, dynamic labels, final eigengenes/kME,
module-trait associations and intramodular connectivity are compared. Correlation/TOM
use labelled row probes: 512,000 elements at 8K and 780,000 at 12K (64 and 65 rows,
respectively). This is not full-square-matrix elementwise validation.

## Benchmark methodology

The pre-fix audit measured a default reference-BLAS environment (A) and a controlled
Accelerate environment (B). A uses R WGCNA's default one thread and Rust's default
available parallelism. Some A repeats overlapped audit workloads and are descriptive.
B configures R WGCNA and BLAS pools to eight threads and Rust faer/Rayon pools to eight;
these are pool settings, not a claim of exactly eight total OS threads.

The non-overlapping timing sum includes input/normalization, MAD, full module detection
and module-trait association. Separately measured cor/TOM are not added again because
module detection recomputes them. Output writing is outside that sum.

| scenario | genes | implementation | repeats | median_seconds | min_seconds | max_seconds | sd_seconds |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | 8000 | r | 3 | 122.16699999999999 | 116.79700000000001 | 173.05100000000004 | 31.04440602320062 |
| A | 8000 | rust | 3 | 16.7989 | 13.559999999999999 | 17.9712 | 2.284857615549235 |
| A | 12000 | r | 1 | 405.866 | 405.866 | 405.866 |  |
| A | 12000 | rust | 1 | 52.1384 | 52.1384 | 52.1384 |  |
| B | 8000 | r | 3 | 9.043 | 9.017 | 9.195 | 0.09614572273377583 |
| B | 8000 | rust | 3 | 10.4537 | 9.793800000000001 | 10.566 | 0.41720741044872756 |
| B | 12000 | r | 3 | 21.389 | 19.796000000000003 | 22.779999999999998 | 1.4931390870690264 |
| B | 12000 | rust | 3 | 34.7317 | 33.932500000000005 | 41.4004 | 4.100402637546706 |


These performance observations precede the correctness fixes and do not benchmark the
new implementation. Three default 8K repeats and one default 12K repeat completed;
two further default 12K repeats were explicitly skipped after the first took about
834 seconds. Controlled B has three repeats at each size. Median/range/standard
deviation and timing boundaries are preserved rather than selecting favourable runs.
These ratios are environment-specific, not algorithm-only speedups.

See [benchmark summaries](../audit/benchmark_summary.csv),
[per-repeat timings](../audit/benchmark_details.csv), and the historical audit's
memory/measurement limitations. Default-process time logs report peak RSS;
B memory sampling is a one-second lower bound, not an exact peak.

## Reproduction and published evidence

Use [validation/README.md](../README.md) for portable commands, dependency requirements,
fresh output directories and the publication manifest. The selected evidence includes
fixtures, numerical summaries, command logs, statistical-grid values and key gene/module
stage exports. Large raw matrices and redundant exports stay local and are catalogued
by checksum; they are not silently represented as published full-matrix evidence.

The source/fixture hashes identify the tested numerical content. Documentation and
publication changes do not alter the tested production algorithms. Published reports
are dated evidence of this run, while the Quarto report reads the result CSVs directly.

## Limits and follow-up coverage

- 8K and 12K are subsets of one biological source, not independent biological datasets.
- The tested configurations do not establish parity for all missing-data patterns,
  block splitting, external branch splitting, all possible powers or every parameter
  interaction. Unsupported or unexecuted capabilities are not counted as passing.
- No post-fix full default-reference-BLAS 12K rerun or performance repetition was done;
  fresh real-data accuracy diagnostics used the stated Accelerate environment.
- Trait input validation is intentionally stricter: ambiguous IDs and rows formerly
  zero-filled or truncated now fail rather than compute silently.
- Complete large adjacency/correlation/TOM square-matrix comparisons remain unexecuted.
  Additional datasets and streamed full-matrix checks are useful next coverage steps.

No scientific thresholds or biological input values were changed to obtain PASS.
