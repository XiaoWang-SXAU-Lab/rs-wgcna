//! Module membership (kME): correlation of each gene with module eigengenes.

use faer::Mat;

use crate::stats::pearson_vec;

/// kME of every gene against each module eigengene.
///
/// `x` is genes x samples; `mes` is samples x modules (each column an
/// eigengene). Returns genes x modules.
pub fn kme(x: &Mat<f64>, mes: &Mat<f64>) -> Mat<f64> {
    let g = x.nrows();
    let n = x.ncols();
    let k = mes.ncols();
    let mut out = Mat::<f64>::zeros(g, k);
    for i in 0..g {
        let gene: Vec<f64> = (0..n).map(|j| x[(i, j)]).collect();
        for mm in 0..k {
            let me: Vec<f64> = (0..n).map(|j| mes[(j, mm)]).collect();
            out[(i, mm)] = pearson_vec(&gene, &me);
        }
    }
    out
}
