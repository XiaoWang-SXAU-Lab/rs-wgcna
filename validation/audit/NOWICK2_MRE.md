# signed Nowick 2 TOM mismatch

Classification: FAIL, not NOT SUPPORTED. Both implementations accept this TOM type. In all six tested combinations (three network types × two denominators), TOM values fail 1e-9 tolerance. Maximum TOM error on the synthetic fixture is about 0.23193. Tree differences follow; in the signed-network cases exact dynamic/final labels are only 20% identical, but ARI is 1 (same partition with different label names). This distinction is explicitly retained in the report.

Minimal R reference:

```r
WGCNA::TOMsimilarityFromExpr(cbind(a=as.numeric(1:4),b=-as.numeric(1:4)), networkType="unsigned", TOMType="signed Nowick 2", power=6, verbose=0)
```

R off-diagonal is +0.5. Rust's public tom_from_cor on correlation [[1,-1],[-1,1]] gives -0.5 with these parameters. Executable reproducer: `cargo test --test audit_numeric_contract audit_signed_nowick2_unsigned_network_matches_reference`. It intentionally fails at the audited commit. The R reference output is `artifacts/nowick2_reference_verified.log`; the full parameter cases are matrix_010/011/022/023/034/035.

Rust `src/tom.rs` unconditionally maps both SignedNowick and SignedNowick2 to UnsignedKeepSign. The observed R `TOMsimilarityFromExpr` type-2 result does not preserve the negative correlation sign in the unsigned network. The exported R wrapper and its accepted TOM types are captured in `artifacts/TOMsimilarityFromExpr_reference.R`. This audit proves the observable discrepancy and identifies the suspect mapping; it does not modify the algorithm or assume all Nowick variants have identical semantics.
