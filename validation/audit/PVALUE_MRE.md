# Student correlation P-value failure

Classification: FAIL. No production algorithm modified.

`cargo test --test audit_numeric_contract audit_near_zero_student_pvalue_matches_reference` is an intentionally failing regression. `examples/audit_stats.rs` exercises the actual public Rust functions on 252 deterministic (r,n) cases, and `stats_reference.R` computes the R WGCNA reference on the same values.

The 252-case Student grid fails in 124 cases at the unchanged absolute tolerance 1e-9; maximum absolute error is 0.999995770203285. At r=1e-8, n=5, R yields approximately 0.9999999872676045 while Rust yields approximately 0.0000042170643195. This can change significance conclusions for weak correlations. Fisher P-values pass this grid (max absolute error 1.7311e-12).

The incomplete beta implementation in `src/stats.rs::betai` always evaluates one continued-fraction branch, with a finite iteration cap and no explicit non-convergence error or complementary transformation near x=1. The numerical failure is proven; identifying a complete stable replacement is outside this audit. Do not label the result as machine precision based on the original five-module synthetic fixture.

Evidence: `artifacts/stats_grid/{all_cases.csv,parity_details.csv,rust.tsv,r.tsv}`. The 8K real-data diagnostic also fails P-value tolerance (max error 4.1016e-8) despite matching modules and correlations; the 12K diagnostic has a much larger P-value discrepancy as well as a module assignment difference.
