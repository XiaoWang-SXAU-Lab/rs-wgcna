<p align="center"><img src="rsWGCNA_icon.png" alt="rsWGCNA logo" width="180"></p>

# rs-WGCNA

**A high-performance, zero-configuration reimplementation of the WGCNA computational kernel in pure Rust.**

[![中文](https://img.shields.io/badge/中文-README.zh--CN-red.svg)](README.zh-CN.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/Rust-1.99%2B-orange.svg)](https://www.rust-lang.org)
[![Parity](https://img.shields.io/badge/Validation-52%20configs%20PASS-brightgreen.svg)](#verification)
[![Zero-config](https://img.shields.io/badge/BLAS-zero--config-lightgrey.svg)](#why-rewrite)
[![linalg](https://img.shields.io/badge/linalg-faer-blueviolet.svg)](https://docs.rs/faer)

> WGCNA computation in pure Rust, with no system BLAS configuration and reproducible comparisons against R WGCNA.

---

## Overview

`rs-wgcna` is a **pure-Rust reimplementation of the computational kernel of the R [`WGCNA`](https://cran.r-project.org/package=WGCNA) package**, covering the full core pipeline of weighted gene co-expression network analysis:

| Step | Module |
|---|---|
| Expression cleaning / filtering (normalization, low-expression filter, MAD/Var gene selection) | `filter` |
| Soft-threshold β selection (scale-free fit) | `sft` |
| Correlation / adjacency matrix | `correlation`, `adjacency` |
| Topological overlap matrix (TOM) | `tom` |
| Hierarchical clustering (UPGMA) + dynamic tree cut | `hclust`, `cutree` |
| Module eigengenes + `mergeCloseModules` | `eigengene`, `merge` |
| Module–trait association + kME/GS | `module_trait`, `kme` |
| End-to-end module detection (single-block module workflow) | `pipeline` |

Linear algebra is powered by [`faer`](https://docs.rs/faer) (pure Rust, SIMD + multithreading), with **no dependence on any system BLAS**.

---

## Why rewrite

1. **No system BLAS setup** — SIMD and parallel linear algebra are provided by `faer`.
2. **Rust integration** — a typed library and worker expose the core single-block workflow.
3. **Explicit resource limits** — current dense network calculations use O(G²) memory.
   Single precision and out-of-core TOM remain [planned features](#roadmap).

---

## How it was reconstructed

The implementation follows the reference WGCNA algorithms. Numerical agreement is assessed with explicit tolerances and identifier checks; tested coverage is listed below.

| Module | Ported from | Key detail |
|---|---|---|
| `tom` | WGCNA C source `src/networkFunctions.c` | Includes the hidden `UnsignedKeepSign` adjacency used when `networkType=unsigned` + `TOMType=signed`; all 6 TOM types |
| `hclust` | R Fortran `src/library/stats/src/hclust.f` (Murtagh algorithm) | Average linkage with reference tie-breaking and tested merge/height agreement |
| `cutree` | `dynamicTreeCut:::cutreeHybrid` (690-line R) | Core-scatter/gap criteria and PAM stage |
| `sft` | `pickSoftThreshold` + `scaleFreeFitIndex` | R 4.5 `cut()` break semantics |
| `eigengene` | `moduleEigengenes` | Standardize + first right singular vector + "along average" sign alignment |
| `merge` | `mergeCloseModules` | relabel + `standardColors` |
| `pipeline` | `blockwiseModules` | kME cleanup + Fisher-p-value gene reassignment (`reassignThreshold=1e-6`) |

---

## Verification

The parity script [`scripts/run_parity.sh`](scripts/run_parity.sh) does everything in one shot: generate data → Rust (release) → R → element-wise comparison (tol = 1e-9).

The default deterministic fixture has 150 genes × 40 samples. Its measured errors are:

| Output | Max diff vs R |
|---|---|
| cor / adjacency / TOM / eigengene / kME / module–trait correlation | ≤5e-15 on the default fixture |
| Soft-threshold fit indices | 9.6e-14 on the default fixture |
| Module–trait p-value | 1.6e-12 on the default fixture |
| Student correlation p-values (252-point grid) | 2.31e-12 |
| Final module colors (real-data diagnostics) | **100%** at both 8k and 12k genes; ARI = 1 |

| Validation | Result |
|---|---|
| Rust tests | 22 PASS |
| Network/module parameter configurations | 52 PASS |
| Comparison-validator controls | 15 PASS |
| Exported comparisons | 664 PASS |

Continuous comparisons use an absolute tolerance of `1e-9`. Gene/sample IDs,
non-finite values and required files are checked. Final labels and ARI cover all
8,000/12,000 selected genes; large correlation/TOM comparisons use labelled row
probes. Runtime and numerical agreement depend on the tested inputs and configuration.

**Scientific audit, fixes and verification: OpenAI Codex (ChatGPT).**
[English audit](validation/fixes/FIX_VERIFICATION.md) ·
[中文审计](validation/fixes/FIX_VERIFICATION.zh-CN.md) ·
[Machine-readable results](validation/README.md).

---

## Benchmark

The included real-data fixture contains 36,953 genes × 56 samples, with MAD-selected
8K/12K subsets. Benchmark results distinguish default reference BLAS from controlled
Accelerate BLAS and record thread settings, timing boundaries and variation. Runtime
ratios describe those environments; no universal or algorithm-only speedup is claimed.

See the [audit benchmark results](validation/audit/benchmark_summary.csv),
[audit methodology](validation/audit/SCIENTIFIC_PARITY_AUDIT.md), and
[data-driven Quarto report](report/benchmark-fixed.html)
([source](report/benchmark.qmd)).

---

## Installation

### As a library (recommended)

Add the git dependency to your `Cargo.toml`:

```toml
[dependencies]
rs-wgcna = { git = "https://github.com/XiaoWang-SXAU-Lab/rs-wgcna", branch = "master" }
```

```rust
use rs_wgcna::io::read_expression;
use rs_wgcna::pipeline::blockwise_modules_single;
use rs_wgcna::network::{NetworkType, TomDenom, TomType};

let exp = read_expression("expr.tsv").unwrap();
let colors = blockwise_modules_single(
    &exp.data, 6.0, NetworkType::Unsigned, TomType::Signed, TomDenom::Min,
    30, 0.25, 2, 0.995, 8,
);
```

### As a binary (parity harness)

```bash
cargo install --git https://github.com/XiaoWang-SXAU-Lab/rs-wgcna
```

> **crates.io**: not published yet — publishing would enable `cargo add rs-wgcna` / `cargo install rs-wgcna` directly (see [Roadmap](#roadmap)).

---

## Usage

### Build

```bash
cargo build --release
```

### Run the parity check

```bash
./scripts/run_parity.sh 6 1e-9
```

### Run the real-data benchmark

```bash
cargo build --release --example real_bench
OUT="$(mktemp -d)"
./target/release/examples/real_bench Demo/WGCNA_fpkm_data.txt Demo/WGCNA_data_trait.txt "$OUT/rust" 8000
Rscript scripts/real_bench.R Demo/WGCNA_fpkm_data.txt Demo/WGCNA_data_trait.txt "$OUT/r" 8000
```

---

## Project structure

```
rs_WGCNA/
├── src/                 # core library
│   ├── correlation.rs   #   Pearson correlation
│   ├── adjacency.rs     #   adjacency
│   ├── tom.rs           #   topological overlap matrix
│   ├── sft.rs           #   soft threshold
│   ├── hclust.rs        #   hierarchical clustering
│   ├── cutree.rs        #   dynamic tree cut
│   ├── eigengene.rs     #   module eigengenes
│   ├── merge.rs         #   mergeCloseModules
│   ├── filter.rs        #   cleaning / filtering
│   ├── module_trait.rs  #   module–trait
│   ├── kme.rs, stats.rs, pipeline.rs, network.rs, io.rs, error.rs
├── examples/            # benchmark & parity entry points
├── scripts/             # R reference & parity scripts
├── Demo/                # real demo data
├── bench_out/           # historical benchmark evidence
├── validation/          # bilingual audits, fixtures, summaries and selected evidence
└── report/              # Quarto report
```

---

## Roadmap

- [x] Core pipeline parity (cor/TOM/hclust/cutree/eigengene/merge/module-trait/kME/GS)
- [x] Real-data benchmark report
- [ ] `f32` single-precision mode (half memory, double SIMD throughput)
- [ ] Tiled / out-of-core TOM (eliminate O(G²) memory, drop `maxBlockSize` entirely)
- [ ] Publish to crates.io
- [ ] Python / R bindings (`pyo3` / `extendr`)

---

## Contributors

- [**Shawn Wang**](https://github.com/ShawnWx2019) — design & implementation
- [**DeepSeek**](https://github.com/deepseek-ai) — AI co-development
- **OpenAI Codex (ChatGPT)** — scientific parity audit, numerical/input fixes, regression tests and bilingual audit reports

## License

[MIT](LICENSE) © 2026 Shawn Wang (Xiao Wang-SXAU-Lab) and contributors
