#!/usr/bin/env bash
# End-to-end parity check: generate data, run Rust and R, compare.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
POWER="${1:-6}"
TOL="${2:-1e-9}"

CARGO="$HOME/.cargo/bin/cargo"

echo "== 1. generate test data =="
Rscript "$HERE/gen_testdata.R" "$ROOT/testdata/expr.tsv"

echo "== 2. build Rust (release) =="
cd "$ROOT"
"$CARGO" build --release --quiet

echo "== 3. run Rust =="
"$ROOT/target/release/rs-wgcna" \
    "$ROOT/testdata/expr.tsv" "$POWER" "$ROOT/rust_out" \
    "$ROOT/testdata/colors.txt" "$ROOT/testdata/trait.tsv"

echo "== 4. run R reference =="
Rscript "$HERE/r_reference.R" \
    "$ROOT/testdata/expr.tsv" "$ROOT/r_out" "$POWER" \
    "$ROOT/testdata/colors.txt" "$ROOT/testdata/trait.tsv"

echo "== 5. compare =="
Rscript "$HERE/compare.R" "$ROOT/rust_out" "$ROOT/r_out" "$TOL"
