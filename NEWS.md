# News / Changelog

All notable changes to rs-WGCNA are documented in this file. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Fixed

- Stabilize Student correlation P-values for weak correlations using the
  complementary incomplete beta calculation rather than an unconverged fraction;
  handle infinite Fisher z statistics at perfect correlation explicitly.
- Preserve fractional soft-threshold powers in adjacency/TOM while retaining the
  integer fast path; match R's network mapping for `signed Nowick 2` TOM.
- Reuse the original dynamic-module eigengenes during post-cleanup reassignment,
  excluding deleted modules and removing modules made too small by reassignment.
- Align uploaded traits by sample ID in the CLI, worker and benchmarks; reject
  duplicate/missing IDs and malformed rows rather than silently pairing values.
- Export benchmark continuous values at full precision for strict comparisons.

### Validation and audit

- A scientific audit of baseline `4f505542` reproduced the 12K gene mismatch and
  identified fractional-power truncation, `signed Nowick 2` mapping differences,
  unstable weak-correlation Student P-values, and positional trait alignment.
- The 12K discrepancy involved `ENSG00000149541`: R retained turquoise while Rust
  assigned yellow. Reassignment used original dynamic-module eigengenes in R and
  recomputed cleaned-module eigengenes in Rust. Matching the eigengene basis
  resolves the difference without changing the `1e-6` reassignment threshold.
- Verification passes: 22 Rust tests, 52 parameter configurations, 15 validator
  controls, 664 exported comparisons, and 252-point Student/Fisher grids at the
  unchanged absolute tolerance `1e-9`. Final labels match R for all 8K/12K genes
  (100%, ARI=1). Large-data correlation/TOM coverage uses labelled row probes.
- Add machine-readable results and bilingual reports:
  [English](validation/fixes/FIX_VERIFICATION.md) /
  [中文](validation/fixes/FIX_VERIFICATION.zh-CN.md).
- Scientific audit, implementation fixes, regression testing and audit reporting
  were completed by **OpenAI Codex (ChatGPT)** at the maintainer's request.
- README validation and performance descriptions now state measured coverage and
  resource requirements; single precision/out-of-core computation remain planned.

## [0.5.0] — 2026-10-26

### Fixed

- Module-detail and hub steps now honor the module-trait "trait source" (the
  uploaded table or the sample_info one-hot column) instead of hardcoding
  `trait.tsv`. GS/MM are now computed against the correct trait.
- Silent fallbacks replaced with loud errors: `trait_column` (no more index-0
  fallback), and network/tom type parsing panic on unknown values instead of
  silently defaulting.

## [0.4.0] — 2026-10-26

### Added

- Module-trait trait source: use the uploaded trait table, or generate a
  one-hot categorical trait matrix from `sample_info.tsv` by a chosen column
  (e.g. `group` / `condition` / `time_h`), aligned to the expression samples.

## [0.3.0] — 2026-10-24

### Changed

- SFT: scan the full 1..20 power range (no early stop) and require a negative
  log-log slope when selecting the soft threshold, avoiding spurious R² passes
  at low power where the network is still dense.
- The network step now honors the user-selected power instead of always using
  the SFT estimate (one-step mode still auto-syncs the SFT power).

### Performance

- Parallelized the per-gene dispersion MLE in `vst` (≈3.4× faster filter on
  large matrices).
- Parallelized the O(g²) TOM adjacency / connectivity / normalization loops.

## [0.2.0] — 2026-10-23

### Added

- **DESeq2 variance-stabilizing transformation (vst), pure Rust** — no R/DESeq2
  dependency. Median-ratio size factors, gene-wise NB dispersion MLE with the
  Cox-Reid correction and the fitDisp revert rule, parametric Gamma-GLM fit,
  and the exact vst formula. Verified against R `DESeq2::vst`: max relative
  difference **~0.09%**, Pearson correlation **0.9999999**.
- `vst` normalization option for raw counts (e.g. RSEM expected counts are
  ceiling'd to integers first).
- CSV input: delimiter auto-detection (tab vs comma) in `io`.
- Extra normalization methods: `log2(x+1)`, `log10(x)`, `log2(x)`.
- Excel workbook input (xlsx/xls/xlsb/ods) via `calamine` (magic-byte sniffing).
- Hub-gene network: hub-and-spoke extraction (hub genes + neighbors, per-hub
  top-N) with per-node kME/GS/degree metrics.

### Changed

- `stats::ln_gamma` is now public (used by the vst module).

## [0.1.0] — 2026-10-21

### Added

- Pure-Rust reimplementation of the WGCNA computational kernel:
  correlation, adjacency, TOM, hierarchical clustering, dynamic tree cut,
  module eigengenes, `mergeCloseModules`, module–trait association, kME/GS,
  and an end-to-end `blockwiseModules`-equivalent pipeline.
- `rs-wgcna-worker` binary (`--job-dir --step`) with steps: `filter`,
  `soft_threshold`, `network`, `module_trait`, `module_detail`, `hub`.
- Worker extras: sample clustering tree, filter report, outlier-sample
  exclusion, experience-power fallback.

### Changed

- Soft-threshold selection ~30× faster at >10k genes (`powi` + rayon +
  parallel connectivity + early stop).
- Filter excludes zero-variance genes to avoid NaN downstream.

### Fixed

- NaN-safe soft-threshold (zero-variance genes no longer panic).

## Verification

Continuous values match R WGCNA to ~1e-14; module colors match **100%** at
8k genes and **99.99%** at 12k genes. See [README](README.md#verification).
