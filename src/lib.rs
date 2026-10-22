//! rs-wgcna: a Rust re-implementation of the WGCNA computational kernel.
//!
//! The goal is numerical parity with the R `WGCNA` package (within tolerance),
//! while removing the R/BLAS configuration burden and memory limits.
//!
//! Convention: an expression matrix is stored as `genes (rows) x samples (cols)`,
//! i.e. `Mat::<f64>` with `nrows() == n_genes` and `ncols() == n_samples`. This
//! matches the on-disk layout used by the reference Shiny app ("gene ID in row,
//! sample ID in column"). WGCNA internally uses samples-x-genes; the translation
//! is handled inside each routine.

pub mod adjacency;
pub mod correlation;
pub mod cutree;
pub mod eigengene;
pub mod error;
pub mod filter;
pub mod hclust;
pub mod io;
pub mod kme;
pub mod merge;
pub mod module_trait;
pub mod network;
pub mod pipeline;
pub mod sft;
pub mod stats;
pub mod tom;
pub mod vst;

pub use error::{Result, WgcnaError};
pub use network::{NetworkType, TomDenom, TomType};
