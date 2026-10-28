# Reproduction and evidence contracts

Run from the repository root. Commands and actual exit codes for this run are in `test_manifest.json`; full stdout/stderr in `test_execution_log.txt` and per-command `artifacts/*.log`. `file_hashes.json` identifies audited production sources, original biological data and old benchmark files. No commit or push was performed.

Prerequisites: Rust/Cargo (lockfile retained), R with WGCNA, dynamicTreeCut, fastcluster, jsonlite, knitr/rmarkdown, Quarto. See `environment.json`. No biological data download or substitution is needed; both inputs come from `Demo/`.

Safe fresh core parity run:

```sh
bash scripts/run_parity.sh 6 1e-9
```

The final printed evidence directory is freshly allocated with mktemp; optional third argument specifies a directory. Do not reuse a valuable prior output directory. This wrapper uses existing gen_testdata, Rust CLI, R reference and strict compare scripts, with unchanged default scientific parameters.

Other executable checks:

```sh
cargo test --workspace
Rscript validation/audit/negative_controls.R .
cargo build --release --example audit_matrix --example audit_real --example audit_stats --example audit_reassignment --example audit_perf
```

`cargo test` is expected to be nonzero at the audited commit: four unignored integration regressions preserve discovered bugs. The original twelve unit tests still pass. Follow the per-case commands in the manifest to reproduce the 52 parameter cases; the same baseline fixture and deterministic challenging fixture are reused. `run_matrix.py` enumerates all cases; it refuses overwriting completed case directories. To run a complete new audit, use a fresh copy of the repository and copy the original fixture-generating scripts from source commit `4f5055425ec9999dfe7b15422d24b5bb72cb1f7e` for the unmodified baseline. Do not relabel tests on the hardened tree as the original baseline.

Optional diagnostic generators (`prepare_matrix.py`, `prepare_real_diagnostics.py`, `prepare_perf.py`) generate TEST-ONLY examples by extending the existing harnesses. They do not alter numerical Rust code. They need not be run to use the checked-in test harness sources. `prepare_matrix.py` intentionally regenerates `scripts/r_reference.R` from the audited HEAD; preserve any subsequent validation edits before regenerating.

The real diagnostic script exports 64 evenly spaced rows of the complete gene correlation/TOM matrix (65 at 12K because of integer stride), full clustering merges/heights, raw dynamic cut, cleanup labels, final eigengenes/kME, connectivity and R premerge labels. Correlation/TOM agreement is therefore sampled, NOT a claim of full 8K/12K element-wise matrix parity. Dynamic labels and final colors are compared for all genes. Unlabelled synthetic colors follow input gene order, proven using each side's labelled correlation matrix. Trait sample IDs must be unique and are explicitly verified against expression IDs.

R real_bench includes the grey eigengene; Rust excludes grey. The comparison explicitly excludes only that documented extra R row/column and records the exclusion. Missing required diagnostic files FAIL. Missing outputs never become PASS. Original benchmarks round Rust module–trait values to six decimals: their 1e-9 FAIL is an output precision limitation; full precision diagnostics separately establish algorithm error. No tolerance has been loosened.

Scenario A: unchanged default R reference BLAS, WGCNA 1 thread; default Rust faer with 8 available threads. Three 8K repeats, one 12K repeat. Some repeats overlapped other audit processes and are descriptive only. Extra default 12K repeats are SKIPPED by `repetition_policy.json`, rather than counted as passed.

Scenario B uses a process-local DYLD_LIBRARY_PATH pointing to `artifacts/accelerate_blas`, whose symlinks resolve to the already installed libRblas.vecLib.dylib. Rscript launches a system shell that loses this override on this host; use the direct R executable as logged. `R_HOME` is set only for that child. `R_environment_accelerate8.json` proves the loaded BLAS and WGCNA thread setting. Both sides are capped/configured to eight threads; timings remain implementation/environment measurements, not algorithm-only speedups. System R libraries and configuration are untouched.

`compare_completed.py` compares existing completed runs and emits per-run CSVs; `aggregate_results.py` assembles the audit CSVs/JSON. `report/benchmark.qmd` reads those generated tables; render with:

```sh
(cd report && quarto render benchmark.qmd --output benchmark-audit.html)
```

The original `report/benchmark.html` and original bench_out evidence remain intact. Additional inference limitations are in SCIENTIFIC_PARITY_AUDIT.md.
