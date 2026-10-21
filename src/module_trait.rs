//! Module-trait relationships and gene significance (GS).

use faer::Mat;

use crate::stats::{cor_pvalue, pearson_vec};

/// Correlation + p-value of each module eigengene (rows = samples) with each
/// trait (rows = samples, aligned to the same sample order).
///
/// Returns (cor, p) matrices, each modules x traits.
pub fn module_trait_cor(me: &Mat<f64>, trait_data: &Mat<f64>) -> (Mat<f64>, Mat<f64>) {
    let ns = me.nrows();
    let nm = me.ncols();
    let nt = trait_data.ncols();
    let mut cor = Mat::<f64>::zeros(nm, nt);
    let mut p = Mat::<f64>::zeros(nm, nt);
    for mi in 0..nm {
        let mev: Vec<f64> = (0..ns).map(|j| me[(j, mi)]).collect();
        for ti in 0..nt {
            let tv: Vec<f64> = (0..ns).map(|j| trait_data[(j, ti)]).collect();
            let r = pearson_vec(&mev, &tv);
            cor[(mi, ti)] = r;
            p[(mi, ti)] = cor_pvalue(r, ns);
        }
    }
    (cor, p)
}

/// Gene significance: correlation of each gene with a single trait.
///
/// `x` is genes x samples; `trait_v` is length = n_samples. Returns length =
/// n_genes.
pub fn gene_significance(x: &Mat<f64>, trait_v: &[f64]) -> Vec<f64> {
    let g = x.nrows();
    let n = x.ncols();
    (0..g)
        .map(|i| {
            let gene: Vec<f64> = (0..n).map(|j| x[(i, j)]).collect();
            pearson_vec(&gene, trait_v)
        })
        .collect()
}
