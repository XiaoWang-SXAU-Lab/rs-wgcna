//! Pearson correlation between genes.

use faer::linalg::matmul::matmul;
use faer::{Accum, Mat, Par};

use crate::error::{Result, WgcnaError};

/// Pearson correlation matrix of the genes (rows) of `x`.
///
/// Equivalent to R's `cor(t(x))` for a samples-x-genes `x`; here `x` is
/// genes-x-samples. Diagonal is exactly 1. Genes with zero variance produce
/// `NaN` correlations (matching R's `cor`, which returns NA).
pub fn pearson(x: &Mat<f64>, threads: usize) -> Result<Mat<f64>> {
    let g = x.nrows();
    let n = x.ncols();
    if n < 2 {
        return Err(WgcnaError::InvalidInput(
            "pearson correlation requires >= 2 samples".into(),
        ));
    }

    // Row means and sample standard deviations.
    let mut means = vec![0.0f64; g];
    let mut sds = vec![0.0f64; g];
    for i in 0..g {
        let mut s = 0.0;
        for j in 0..n {
            s += x[(i, j)];
        }
        let m = s / n as f64;
        means[i] = m;
        let mut ss = 0.0;
        for j in 0..n {
            let d = x[(i, j)] - m;
            ss += d * d;
        }
        sds[i] = (ss / (n as f64 - 1.0)).sqrt();
    }

    // Standardize rows.
    let mut z = Mat::<f64>::zeros(g, n);
    for i in 0..g {
        let m = means[i];
        let sd = sds[i];
        for j in 0..n {
            z[(i, j)] = if sd == 0.0 {
                f64::NAN
            } else {
                (x[(i, j)] - m) / sd
            };
        }
    }

    // cor = Z Z^T / (n - 1).
    let mut cor = Mat::<f64>::zeros(g, g);
    let zt = z.transpose();
    let par = Par::rayon(threads);
    matmul(
        &mut cor,
        Accum::Replace,
        &z,
        zt,
        1.0 / (n as f64 - 1.0),
        par,
    );

    // Fix diagonal (and avoid NaN on diagonal for constant genes, matching the
    // fact that cor of a gene with itself is 1 whenever variance > 0).
    for i in 0..g {
        cor[(i, i)] = 1.0;
    }

    Ok(cor)
}
