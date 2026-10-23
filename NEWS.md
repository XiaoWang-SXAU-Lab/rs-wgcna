# News / Changelog

All notable changes to rs-WGCNA are documented in this file. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

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
