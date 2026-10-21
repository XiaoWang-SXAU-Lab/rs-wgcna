# rs-WGCNA

**用纯 Rust 重写的 WGCNA 计算内核——高性能、零 BLAS 配置。**

[![English](https://img.shields.io/badge/English-README-blue.svg)](README.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/Rust-1.99%2B-orange.svg)](https://www.rust-lang.org)
[![Parity](https://img.shields.io/badge/Parity-R%20WGCNA-~100%25-brightgreen.svg)](#验证)
[![Zero-config](https://img.shields.io/badge/BLAS-zero--config-lightgrey.svg)](#为什么重写)
[![linalg](https://img.shields.io/badge/linalg-faer-blueviolet.svg)](https://docs.rs/faer)

> 让 WGCNA 摆脱 R 的性能与内存瓶颈——零 BLAS 配置、结果与 R 版数值等价。

---

## 概述

`rs-wgcna` 是 **R 版 [`WGCNA`](https://cran.r-project.org/package=WGCNA) 包计算内核的纯 Rust 重实现**，覆盖加权基因共表达网络分析的完整核心管线：

| 步骤 | 模块 |
|---|---|
| 表达矩阵清洗 / 过滤（归一化、低表达过滤、MAD/Var 选基因） | `filter` |
| 软阈值 β 选择（scale-free 拟合） | `sft` |
| 相关矩阵 / 邻接矩阵 | `correlation`, `adjacency` |
| 拓扑重叠矩阵 TOM | `tom` |
| 层次聚类（UPGMA）+ 动态切树 | `hclust`, `cutree` |
| 模块 eigengene + `mergeCloseModules` | `eigengene`, `merge` |
| 模块-性状关联 + kME/GS | `module_trait`, `kme` |
| 端到端模块检测（`blockwiseModules` 单 block 等价） | `pipeline` |

线性代数使用 [`faer`](https://docs.rs/faer)（纯 Rust，SIMD + 多线程），**不依赖任何系统 BLAS**。

---

## 为什么重写

解决 R 版 WGCNA 的三个核心痛点：

1. **矩阵计算慢** — R 发行版默认使用 reference-BLAS（单线程、无 SIMD）。普通用户不会去换 OpenBLAS/MKL，导致 `cor` 与 `TOM` 极其缓慢。
2. **零配置** — `faer` 是纯 Rust 实现，优化**内置在二进制里**，用户无需安装/配置任何 BLAS。
3. **内存限制** — R 需要用户手动指定 `maxBlockSize` 分块来绕过 O(G²) 内存。rs-WGCNA 通过 f32 与分块/外存计算消除这一限制（见 [Roadmap](#roadmap)）。

---

## 如何重构的

核心原则是 **“忠实移植”（faithful porting）**——不是重新发明算法，而是逐行翻译 WGCNA 的原始源码，从而保证数值上逐位等价。

| 模块 | 移植来源 | 关键细节 |
|---|---|---|
| `tom` | WGCNA C 源码 `src/networkFunctions.c` | 含默认路径 `networkType=unsigned` + `TOMType=signed` 时内部使用 **`UnsignedKeepSign`** 邻接的隐藏行为；支持全部 6 种 TOM 类型 |
| `hclust` | R Fortran 源码 `src/library/stats/src/hclust.f`（Murtagh 算法） | tie-breaking 用严格 `<`/`>`；`merge`/`height`/`order` 与 R 逐位一致 |
| `cutree` | `dynamicTreeCut:::cutreeHybrid`（690 行 R） | core-scatter/gap 判据与 PAM 阶段 |
| `sft` | `pickSoftThreshold` + `scaleFreeFitIndex` | 复刻 R 4.5 `cut()` 的断点语义 |
| `eigengene` | `moduleEigengenes` | 标准化 + SVD 第一右奇异向量 + “along average” 符号对齐 |
| `merge` | `mergeCloseModules` | relabel + `standardColors` 标准配色 |
| `pipeline` | `blockwiseModules` | kME 清洗 + **Fisher p-value 基因重分配**（`reassignThreshold=1e-6`） |

---

## 验证

对拍脚本 [`scripts/run_parity.sh`](scripts/run_parity.sh) 一键完成：生成合成数据 → Rust（release）→ R → 逐项比较（tol = 1e-9）。

| 输出 | 与 R 的最大差 |
|---|---|
| cor / adjacency / TOM / eigengene / kME / 模块-性状相关 | ~1e-15（机器精度） |
| 软阈值 fitIndices | 9.6e-14 |
| 模块-性状 p 值 | 1.6e-12 |
| 最终模块颜色（离散标签） | 8k 基因 **100%**，12k 基因 **99.99%**（1 个边界基因） |

> 数值量达到机器精度；离散的模块标签在更大规模下可能命中个别恰好落在浮点边界的基因（属固有现象，非逻辑错误）。

---

## 基准

真实乳腺癌 RNA-seq 数据（FPKM，36,953 基因 × 56 样本，5 个分子亚型，见 [`Demo/`](Demo/)）。环境：**Apple M2，8 核**；R 用默认 reference-BLAS。

完整 Quarto 报告见 [`report/benchmark.html`](report/benchmark.html)（源码 [`report/benchmark.qmd`](report/benchmark.qmd)）。

### 8,000 基因（MAD top 8000）

| 步骤 | R（秒） | rs-WGCNA（秒） | 加速 |
|---|---|---|---|
| cor | 1.20 | 0.084 | **14×** |
| TOM | 110.0 | 6.76 | **16×** |
| 模块检测（总） | 116.7 | 11.5 | **10×** |
| **端到端** | **~119** | **~11.6** | **~10×** |

### 12,000 基因（MAD top 12000）

| 步骤 | R（秒） | rs-WGCNA（秒） | 加速 |
|---|---|---|---|
| cor | 2.59 | 0.17 | **15×** |
| TOM | 398.6 | 21.0 | **19×** |
| 模块检测（总） | 406.7 | 35.4 | **11.5×** |
| **端到端** | **~410** | **~35.5** | **~11.5×** |

---

## 安装

### 作为库使用（推荐）

在 `Cargo.toml` 中添加 git 依赖：

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

### 作为二进制（对拍工具）

```bash
cargo install --git https://github.com/XiaoWang-SXAU-Lab/rs-wgcna
```

> **crates.io**：尚未发布——发布后即可直接 `cargo add rs-wgcna` / `cargo install rs-wgcna`（见 [Roadmap](#roadmap)）。

---

## 使用说明

### 构建

```bash
cargo build --release
```

### 跑对拍

```bash
./scripts/run_parity.sh 6 1e-9
```

### 跑真实数据基准

```bash
cargo build --release --example real_bench
./target/release/examples/real_bench Demo/WGCNA_fpkm_data.txt Demo/WGCNA_data_trait.txt bench_out/rust 8000
Rscript scripts/real_bench.R Demo/WGCNA_fpkm_data.txt Demo/WGCNA_data_trait.txt bench_out/r 8000
```

---

## 项目结构

```
rs_WGCNA/
├── src/                 # 计算内核（12 个模块）
│   ├── correlation.rs   #   Pearson 相关
│   ├── adjacency.rs     #   邻接矩阵
│   ├── tom.rs           #   拓扑重叠矩阵
│   ├── sft.rs           #   软阈值
│   ├── hclust.rs        #   层次聚类
│   ├── cutree.rs        #   动态切树
│   ├── eigengene.rs     #   模块 eigengene
│   ├── merge.rs         #   mergeCloseModules
│   ├── filter.rs        #   清洗/过滤
│   ├── module_trait.rs  #   模块-性状
│   ├── kme.rs, stats.rs, pipeline.rs, network.rs, io.rs, error.rs
├── examples/            # 基准/对拍入口
├── scripts/             # R 参考与对拍脚本
├── Demo/                # 真实演示数据
├── bench_out/           # 基准结果（供报告使用）
└── report/              # Quarto 报告
```

---

## Roadmap

- [x] 核心管线对拍通过（cor/TOM/hclust/cutree/eigengene/merge/module-trait/kME/GS）
- [x] 真实数据基准报告
- [ ] `f32` 单精度模式（内存减半、吞吐翻倍）
- [ ] 分块/外存 TOM（消除 O(G²) 内存，彻底去掉 `maxBlockSize`）
- [ ] 发布到 crates.io
- [ ] Python / R 绑定（`pyo3` / `extendr`）

---

## 贡献者

- [**Shawn Wang**](https://github.com/ShawnWx2019) — 设计与实现
- [**DeepSeek**](https://github.com/deepseek-ai) — AI 协同开发

## License

[MIT](LICENSE) © 2026 Shawn Wang (Xiao Wang-SXAU-Lab) and contributors
