# 12K mismatch: evidence against a simple floating-point explanation

Gene: ENSG00000149541. R final turquoise; Rust final yellow. Dynamic tree cut label 1 in both; cleanup label 1 in both; R premerge turquoise.

Source inspection: R `blockwiseModules` retains `AllMEs` calculated from initial dynamic modules when performing the global Fisher-ratio reassignment check. Rust `reassign_by_kme` in `src/pipeline.rs` recomputes eigengenes from the labels AFTER cleanup. These are different scientific inputs to the reassignment test. The R function body captured from installed WGCNA is in `artifacts/blockwiseModules_reference.R`.

Reference recomputation on identical exported labels shows:

| Eigengene basis | Own | Best | P(best)/P(own) | Ratio minus 1e-6 | Decision |
|---|---:|---:|---:|---:|---|
| original dynamic modules | 0.5124264102 | 0.7201410897 | 1.0180061545e-6 | +1.80061545e-8 | retain |
| cleaned modules | 0.5124214676 | 0.7211448429 | 9.1835531036e-7 | -8.16446896e-8 | reassign |

The best kME changes by about 0.001, far beyond roundoff. The ratios lie on opposite sides of the production 1e-6 threshold. Raw-cut and cleaned labels match between implementations; the R actual blockwise cut also matches the independent cut export. Thus the earliest evidenced divergence is the reassignment decision after cleanup, consistent with an algorithmic difference in eigengene basis, not merely color renaming, eigengene sign, tree ties, or machine precision. Direct Rust kME/Fisher measurements on these saved stages are exported separately in `reassignment_rust_measurements_verified.tsv`.

Caveat: Rust's private reassignment function exposes no premerge state. Its state was not instrumented or modified. The conclusion combines actual stage outputs, direct measurements through exported numerical APIs, and source-level dataflow. Independent premerge tracing would make the causal chain even stronger.

All final-module kMEs and relevant stage labels for the discordant gene are in `artifacts/diagnostic_12000_comparison/mismatch_investigation.tsv`; threshold evidence is in `reassignment_thresholds.tsv` and `reassignment_rust_measurements_verified.tsv`. The diagnostic R process uses verified Accelerate BLAS with a one-thread cap; default-reference-BLAS repetitions are separately recorded. Original historical evidence remains unchanged.
