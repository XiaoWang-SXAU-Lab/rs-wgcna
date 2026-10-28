# Fractional power discrepancy

Classification: FAIL for the accepted `f64` parameter 2.5; no production changes made.

Reproduce from repository root:

```sh
cargo build --release --example audit_matrix
B=validation/audit/artifacts/baseline/testdata
D=$(mktemp -d)
target/release/examples/audit_matrix "$B/expr.tsv" 2.5 "$D/rust" "$B/colors.txt" "$B/trait.tsv" unsigned signed min 2 0.15
Rscript scripts/r_reference.R "$B/expr.tsv" "$D/r" 2.5 "$B/colors.txt" "$B/trait.tsv" unsigned signed min 2 0.15
Rscript scripts/compare.R "$D/rust" "$D/r" 1e-9 "$D/comparison" "$B/expr.tsv"
```

`src/adjacency.rs` and `src/tom.rs::adj_transform` cast power to i32 before `powi`, without rejecting fractional inputs. R applies the actual real-valued exponent. Both public Rust functions and the existing harness accept an f64. Integer-only intent is described in a source comment but not enforced. This is an interface/algorithm-contract discrepancy, not floating-point rounding.

Full generated evidence: `artifacts/matrix_038/comparison/parity_details.csv` (adjacency max absolute error 0.0704823914619734, TOM 0.08142744815293). Integer powers 1, 6, 12 pass on this fixture. SFT uses a separate integer power vector, so its PASS does not validate fractional exponent handling.

An independent two-gene example: unsigned correlation 0.25, power 2.5 gives Rust 0.0625 (power 2) versus R 0.03125. Follow-up should either implement real powers or reject unsupported inputs explicitly; this audit does neither.
