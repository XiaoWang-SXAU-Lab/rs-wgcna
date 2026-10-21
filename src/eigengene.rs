//! Module eigengenes.
//!
//! Faithful port of WGCNA's `moduleEigengenes`: standardize each module gene to
//! mean 0 / sd 1, take the first right singular vector of the standardized
//! module matrix (via faer's SVD, matching LAPACK), then align its sign with
//! the average expression ("along average").

use faer::linalg::solvers::Svd;
use faer::Mat;

use crate::stats::pearson_vec;

/// Eigengene (length = n_samples) for the genes in `module`.
///
/// `module` holds row indices into `x` (genes x samples).
pub fn module_eigengene(x: &Mat<f64>, module: &[usize]) -> Vec<f64> {
    let n = x.ncols();
    if module.is_empty() {
        return vec![0.0; n];
    }
    let m = module.len();

    // Standardize each module gene; accumulate per-sample mean of the
    // standardized genes (the "average expression" used for sign alignment).
    let mut z = Mat::<f64>::zeros(m, n);
    let mut aver = vec![0.0f64; n];
    for (r, &gi) in module.iter().enumerate() {
        let mut mean = 0.0;
        for j in 0..n {
            mean += x[(gi, j)];
        }
        mean /= n as f64;
        let mut ss = 0.0;
        for j in 0..n {
            let d = x[(gi, j)] - mean;
            ss += d * d;
        }
        let sd = (ss / (n as f64 - 1.0)).sqrt();
        for j in 0..n {
            let v = if sd == 0.0 {
                0.0
            } else {
                (x[(gi, j)] - mean) / sd
            };
            z[(r, j)] = v;
            aver[j] += v;
        }
    }
    for j in 0..n {
        aver[j] /= m as f64;
    }

    // First right singular vector via thin SVD (V is n x min(m,n)).
    let svd = match Svd::new_thin(z.as_ref()) {
        Ok(s) => s,
        Err(_) => return vec![0.0; n],
    };
    let v = svd.V();
    let mut pc: Vec<f64> = (0..n).map(|j| v[(j, 0)]).collect();

    // Sign alignment "along average".
    let cor_ave = pearson_vec(&aver, &pc);
    if cor_ave.is_finite() && cor_ave < 0.0 {
        for v in pc.iter_mut() {
            *v = -*v;
        }
    }
    pc
}
