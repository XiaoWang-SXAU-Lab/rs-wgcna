# rs-WGCNA

**用纯 Rust 重写的 WGCNA 计算内核——高性能、零 BLAS 配置。**

[![English](https://img.shields.io/badge/English-README-blue.svg)](README.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/Rust-1.99%2B-orange.svg)](https://www.rust-lang.org)
[![Parity](https://img.shields.io/badge/Validation-52%20configs%20PASS-brightgreen.svg)](#验证)
[![Zero-config](https://img.shields.io/badge/BLAS-zero--config-lightgrey.svg)](#为什么重写)
[![linalg](https://img.shields.io/badge/linalg-faer-blueviolet.svg)](https://docs.rs/faer)

> 纯 Rust 的 WGCNA 计算内核，无需配置系统 BLAS，提供与 R WGCNA 的可复现对照验证。

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
| 端到端模块检测（单 block 模块检测） | `pipeline` |

线性代数使用 [`faer`](https://docs.rs/faer)（纯 Rust，SIMD + 多线程），**不依赖任何系统 BLAS**。

---

## 为什么重写

1. **无需配置系统 BLAS** — `faer` 提供 SIMD 和并行线性代数。
2. **Rust 集成** — 类型化库与 worker 提供核心单 block 工作流。
3. **明确资源需求** — 当前稠密网络计算需要 O(G²) 内存；单精度和外存 TOM
   仍在 [Roadmap](#roadmap) 中。

---

## 如何重构的

实现遵循参考 WGCNA 算法，通过明确的数值容差和标识符检查验证一致性；已测试范围见下文。

| 模块 | 移植来源 | 关键细节 |
|---|---|---|
| `tom` | WGCNA C 源码 `src/networkFunctions.c` | 含默认路径 `networkType=unsigned` + `TOMType=signed` 时内部使用 **`UnsignedKeepSign`** 邻接的隐藏行为；支持全部 6 种 TOM 类型 |
| `hclust` | R Fortran 源码 `src/library/stats/src/hclust.f`（Murtagh 算法） | 平均连接、参考 tie-breaking；验证 merge/height 一致性 |
| `cutree` | `dynamicTreeCut:::cutreeHybrid`（690 行 R） | core-scatter/gap 判据与 PAM 阶段 |
| `sft` | `pickSoftThreshold` + `scaleFreeFitIndex` | 复刻 R 4.5 `cut()` 的断点语义 |
| `eigengene` | `moduleEigengenes` | 标准化 + SVD 第一右奇异向量 + “along average” 符号对齐 |
| `merge` | `mergeCloseModules` | relabel + `standardColors` 标准配色 |
| `pipeline` | `blockwiseModules` | kME 清洗 + **Fisher p-value 基因重分配**（`reassignThreshold=1e-6`） |

---

## 验证

对拍脚本 [`scripts/run_parity.sh`](scripts/run_parity.sh) 一键完成：生成合成数据 → Rust（release）→ R → 逐项比较（tol = 1e-9）。

默认确定性 fixture 为 150 基因 × 40 样本，实测误差如下：

| 输出 | 与 R 的最大差 |
|---|---|
| cor / adjacency / TOM / eigengene / kME / 模块-性状相关 | 默认 fixture ≤5e-15 |
| 软阈值 fitIndices | 默认 fixture 9.6e-14 |
| 模块-性状 p 值 | 默认 fixture 1.6e-12 |
| Student 相关 p 值（252 点网格） | 2.31e-12 |
| 最终模块颜色（真实数据诊断） | 8k 和 12k 基因均 **100%**，ARI = 1 |

| 验证项目 | 结果 |
|---|---|
| Rust 测试 | 22 PASS |
| 网络/模块参数配置 | 52 PASS |
| 比较验证器控制测试 | 15 PASS |
| 已导出结果比较 | 664 PASS |

连续量采用绝对容差 `1e-9`，检查基因/样本 ID、非有限值和必需文件。
8K/12K 的最终标签与 ARI 覆盖全部保留基因；大数据 cor/TOM 使用带 ID 的行探针。
耗时和数值一致性以已测试的输入与配置为范围。

**科学审计、修复和复测由 OpenAI Codex（ChatGPT）完成。**
[中文审计](validation/fixes/FIX_VERIFICATION.zh-CN.md) ·
[English audit](validation/fixes/FIX_VERIFICATION.md) ·
[机器可读结果](validation/README.md)。

---

## 基准

仓库内真实数据 fixture 包含 36,953 基因 × 56 样本，使用 MAD 选择 8K/12K 子集。
基准结果区分默认 reference BLAS 和受控 Accelerate BLAS，记录线程设置、计时边界和
重复测量波动。耗时比值仅描述相应环境，不宣称普遍或纯算法加速。

详见[审计基准结果](validation/audit/benchmark_summary.csv)、
[审计方法](validation/audit/SCIENTIFIC_PARITY_AUDIT.md)及
[数据驱动 Quarto 报告](report/benchmark-fixed.html)（[源码](report/benchmark.qmd)）。

---

## 安装

### 作为库使用（推荐）

在 `Cargo.toml` 中添加 git 依赖：

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
OUT="$(mktemp -d)"
./target/release/examples/real_bench Demo/WGCNA_fpkm_data.txt Demo/WGCNA_data_trait.txt "$OUT/rust" 8000
Rscript scripts/real_bench.R Demo/WGCNA_fpkm_data.txt Demo/WGCNA_data_trait.txt "$OUT/r" 8000
```

---

## 项目结构

```
rs_WGCNA/
├── src/                 # 计算内核
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
├── bench_out/           # 历史基准证据
├── validation/          # 中英文审计、fixture、汇总与选定证据
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
- **OpenAI Codex（ChatGPT）** — 科学一致性审计、数值与输入修复、回归测试和中英文审计报告

## License

[MIT](LICENSE) © 2026 Shawn Wang (Xiao Wang-SXAU-Lab) and contributors
