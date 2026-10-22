<p align="center"><img src="rsWGCNA_icon.png" alt="rsWGCNA logo" width="180"></p>

# rs-WGCNA

**A high-performance, zero-configuration reimplementation of the WGCNA computational kernel in pure Rust.**

[![中文](https://img.shields.io/badge/中文-README.zh--CN-red.svg)](README.zh-CN.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/Rust-1.99%2B-orange.svg)](https://www.rust-lang.org)
[![Parity](https://img.shields.io/badge/Parity-R%20WGCNA-~100%25-brightgreen.svg)](#verification)
[![Zero-config](https://img.shields.io/badge/BLAS-zero--config-lightgrey.svg)](#why-rewrite)
[![linalg](https://img.shields.io/badge/linalg-faer-blueviolet.svg)](https://docs.rs/faer)

> Free WGCNA from R's performance and memory bottlenecks — zero BLAS configuration, numerically equivalent to the R package.

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
| End-to-end module detection (single-block `blockwiseModules` equivalent) | `pipeline` |

Linear algebra is powered by [`faer`](https://docs.rs/faer) (pure Rust, SIMD + multithreading), with **no dependence on any system BLAS**.

---

## Why rewrite

It addresses three core pain points of the R WGCNA:

1. **Slow matrix ops** — R ships with the reference BLAS (single-threaded, no SIMD). Most users never switch to OpenBLAS/MKL, making `cor` and `TOM` extremely slow.
2. **Zero configuration** — `faer` is pure Rust; the optimization is baked into the binary. Users install/configure nothing.
3. **Memory limits** — R forces users to hand-tune `maxBlockSize` to work around O(G²) memory. rs-WGCNA removes this via f32 and out-of-core computation (see [Roadmap](#roadmap)).

---

## How it was reconstructed

The guiding principle is **faithful porting** — not reinventing the algorithms, but translating WGCNA's original source line-by-line to guarantee bit-level numerical equivalence.

| Module | Ported from | Key detail |
|---|---|---|
| `tom` | WGCNA C source `src/networkFunctions.c` | Includes the hidden `UnsignedKeepSign` adjacency used when `networkType=unsigned` + `TOMType=signed`; all 6 TOM types |
| `hclust` | R Fortran `src/library/stats/src/hclust.f` (Murtagh algorithm) | Strict `<`/`>` tie-breaking; `merge`/`height`/`order` identical to R |
| `cutree` | `dynamicTreeCut:::cutreeHybrid` (690-line R) | Core-scatter/gap criteria and PAM stage |
| `sft` | `pickSoftThreshold` + `scaleFreeFitIndex` | R 4.5 `cut()` break semantics |
| `eigengene` | `moduleEigengenes` | Standardize + first right singular vector + "along average" sign alignment |
| `merge` | `mergeCloseModules` | relabel + `standardColors` |
| `pipeline` | `blockwiseModules` | kME cleanup + Fisher-p-value gene reassignment (`reassignThreshold=1e-6`) |

---

## Verification

The parity script [`scripts/run_parity.sh`](scripts/run_parity.sh) does everything in one shot: generate data → Rust (release) → R → element-wise comparison (tol = 1e-9).

| Output | Max diff vs R |
|---|---|
| cor / adjacency / TOM / eigengene / kME / module–trait correlation | ~1e-15 (machine precision) |
| Soft-threshold fit indices | 9.6e-14 |
| Module–trait p-value | 1.6e-12 |
| Final module colors (discrete labels) | **100%** at 8k genes, **99.99%** at 12k genes (1 boundary gene) |

> Continuous values reach machine precision; discrete module labels may hit an occasional gene sitting exactly on a floating-point boundary at larger scale (inherent, not a logic bug).

---

## Benchmark

Real breast-cancer RNA-seq data (FPKM, 36,953 genes × 56 samples, 5 molecular subtypes; see [`Demo/`](Demo/)). Environment: **Apple M2, 8 cores**; R uses the default reference-BLAS.

The full Quarto report is at [`report/benchmark.html`](report/benchmark.html) (source: [`report/benchmark.qmd`](report/benchmark.qmd)).

### 8,000 genes (MAD top 8000)

| Step | R (s) | rs-WGCNA (s) | Speedup |
|---|---|---|---|
| cor | 1.20 | 0.084 | **14×** |
| TOM | 110.0 | 6.76 | **16×** |
| module detection (total) | 116.7 | 11.5 | **10×** |
| **end-to-end** | **~119** | **~11.6** | **~10×** |

### 12,000 genes (MAD top 12000)

| Step | R (s) | rs-WGCNA (s) | Speedup |
|---|---|---|---|
| cor | 2.59 | 0.17 | **15×** |
| TOM | 398.6 | 21.0 | **19×** |
| module detection (total) | 406.7 | 35.4 | **11.5×** |
| **end-to-end** | **~410** | **~35.5** | **~11.5×** |

---

## Installation

### As a library (recommended)

Add the git dependency to your `Cargo.toml`:

```toml
[dependencies]
rs-wgcna = { git = "https://github.com/XiaoWang-SXAU-Lab/rs-wgcna", tag = "v0.1.0" }
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
./target/release/examples/real_bench Demo/WGCNA_fpkm_data.txt Demo/WGCNA_data_trait.txt bench_out/rust 8000
Rscript scripts/real_bench.R Demo/WGCNA_fpkm_data.txt Demo/WGCNA_data_trait.txt bench_out/r 8000
```

---

## Project structure

```
rs_WGCNA/
├── src/                 # core library (12 modules)
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
├── bench_out/           # benchmark outputs (used by the report)
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

## License

[MIT](LICENSE) © 2026 Shawn Wang (Xiao Wang-SXAU-Lab) and contributors
