# Trait sample ordering discrepancy

The production parity harness `src/main.rs` and real benchmark `examples/real_bench.rs` read trait sample IDs but discard them. Values are paired positionally with expression samples. R references reorder rows by sample ID. Correct results therefore depend on an undocumented positional input constraint.

Reproduction commands and exit codes are stored under `trait_reordered_*` in `test_manifest.json`; the fixture is the unchanged seed-42 expression matrix and the same trait values with rows reversed. Both inputs retain correct unique sample IDs. Inspect `artifacts/trait_reordered/comparison/parity_details.csv`.

This is a sample alignment defect in the Rust harness/integration, not evidence that the underlying `module_trait_cor` numerical function is wrong. The baseline biological fixture is explicitly verified to have identical sample ordering (`artifacts/input_alignment.log`), so its agreement is not explained by this defect. This audit does not change production input handling. The diagnostic test harness asserts the positional ordering rather than silently assuming it.
