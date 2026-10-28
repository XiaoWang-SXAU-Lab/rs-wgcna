#!/usr/bin/env bash
# End-to-end parity check: generate data, run Rust and R, compare.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
WORK="${3:-$(mktemp -d "${TMPDIR:-/tmp}/rs-wgcna-parity.XXXXXX")}"
mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
echo "Evidence directory: $WORK"
POWER="${1:-6}"
TOL="${2:-1e-9}"

CARGO="$HOME/.cargo/bin/cargo"

echo "== 1. generate test data =="
Rscript "$HERE/gen_testdata.R" "$WORK/testdata/expr.tsv"

echo "== 2. build Rust (release) =="
cd "$ROOT"
"$CARGO" build --release --quiet

echo "== 3. run Rust =="
"$ROOT/target/release/rs-wgcna" \
    "$WORK/testdata/expr.tsv" "$POWER" "$WORK/rust_out" \
    "$WORK/testdata/colors.txt" "$WORK/testdata/trait.tsv"

echo "== 4. run R reference =="
Rscript "$HERE/r_reference.R" \
    "$WORK/testdata/expr.tsv" "$WORK/r_out" "$POWER" \
    "$WORK/testdata/colors.txt" "$WORK/testdata/trait.tsv"

echo "== 5. compare =="
Rscript "$HERE/compare.R" "$WORK/rust_out" "$WORK/r_out" "$TOL" "$WORK/comparison" "$WORK/testdata/expr.tsv"
