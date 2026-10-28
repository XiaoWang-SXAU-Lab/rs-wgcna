# Scientific validation evidence

[Latest English audit](fixes/FIX_VERIFICATION.md) · [最新中文审计](fixes/FIX_VERIFICATION.zh-CN.md)

Scientific auditing, implementation fixes, regression tests and audit reporting were
performed by **OpenAI Codex (ChatGPT)** on 2026-10-28, at the maintainer's request.

## Results

The post-fix run reports 22 Rust tests, 52 parameter configurations, 15 validator
controls and 664 exported comparisons passing at absolute tolerance `1e-9`.
Both 8K/12K final partitions have exact colour agreement of 100% and ARI=1.
Large correlation/TOM comparisons use labelled row probes, not full square matrices.
The original pre-fix FAIL assessment and its failed cases remain separate in `audit/`.

| Location | Contents |
|---|---|
| `fixes/FIX_VERIFICATION.md` / `.zh-CN.md` | Latest English / Chinese scientific audit and fix verification |
| `fixes/test_counts.csv`, `parameter_coverage.csv`, `parity_details.csv`, `module_summary.csv` | Current machine-readable results |
| `fixes/test_manifest.json`, `source_hashes.json` | Exact executed commands, outcomes and tested source/input hashes |
| `audit/SCIENTIFIC_PARITY_AUDIT.md` | Original Chinese pre-fix audit, findings and limitations |
| `audit/benchmark_summary.csv`, `benchmark_details.csv` | Environment-specific pre-fix timings and repetition statistics |
| `audit/artifacts/`, `fixes/artifacts/` | Selected small review artifacts and key module-stage exports |
| `PUBLISHED_EVIDENCE.json` | Checksums, sizes and selection status for retained artifacts |
| `PUBLICATION_CHECKS.json` | Clean-export smoke tests and final publication checks |

## Reproduction

Run from the repository root. Rust/Cargo and R with WGCNA 1.74, dynamicTreeCut,
fastcluster, jsonlite, knitr/rmarkdown are needed; Quarto is needed only for rendering.
Keep the supplied Cargo lockfile. The archived run used Rust 1.99.1 and R 4.5.1.

```sh
cargo test --workspace
bash scripts/run_parity.sh 6 1e-9
```

The parity wrapper creates a fresh temporary output directory. For all 52 configurations
and the real-data diagnostics, choose a new, nonexistent output directory:

```sh
python3 validation/fixes/run_validation.py --output-dir /tmp/rs-wgcna-verification
python3 validation/fixes/summarize.py --output-dir /tmp/rs-wgcna-verification
```

The runner uses the installed `Rscript` and its configured BLAS. It does not depend on
local absolute symlinks from the archived Mac environment. A default reference BLAS can
make the 12K run take many minutes. BLAS/thread changes can affect low-order floating
point results and must be recorded when comparing a new run to the archived one.
The dated audit describes the exact process-local Accelerate environment actually used.

```sh
(cd report && quarto render benchmark.qmd --output benchmark-fixed.html)
```

Quarto reads the published CSVs directly. Summarization writes a separate generated
Chinese draft; it does not overwrite the reviewed bilingual audit reports.
To reproduce the pre-fix failures, use source commit
`4f5055425ec9999dfe7b15422d24b5bb72cb1f7e` and the archived test/reference artifacts.
The historical audit runner is not a claim that the current fixed tree is that baseline.

## Publication policy

Published evidence includes:

- Audit reports, source/command manifests, fixtures and comparison summaries.
- Original and corrected statistical grids, controls and failure reproductions.
- Key real-data gene IDs, module assignments, eigengenes/kME, tree, dynamic/cleanup
  states, connectivity and per-step timings.
- Command output logs and the installed R reference wrapper snapshots.

Large raw correlation/TOM probes, redundant per-case matrix exports, intermediate
render failures and host-specific symlinks remain local. They were not deleted.
`PUBLISHED_EVIDENCE.json` records file hashes and whether an artifact was selected for
Git. Re-running the recorded computations is required to recompute errors from raw
matrices not included in Git; a summary CSV alone is not full-matrix evidence.
The tracked data do not include build targets, credentials or local runtime symlinks.

## 中文说明

最新审计提供中英文版本，并注明由 OpenAI Codex（ChatGPT）完成。
修复前历史证据与修复后结果分开保存；缺失、失败、跳过和未执行不计为通过。
Git 提供汇总、命令、fixture、关键阶段数据与校验值；大矩阵、冗余导出和本机符号链接
保留在本地，不作为已发布的完整平方矩阵证据。完整复测使用上面的 `--output-dir` 指定
新目录，不覆盖已归档结果。当前 runner 使用本机 Rscript 的 BLAS；复现时须记录环境。
