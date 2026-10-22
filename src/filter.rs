//! Expression matrix cleaning / filtering (Step 1 of the reference pipeline):
//! normalization, low-expression filtering, and top-N gene selection by MAD or
//! variance. These are pre-processing steps applied before the WGCNA core.

use faer::Mat;

use crate::error::{Result, WgcnaError};

/// Normalization methods supported by the reference app.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalizeMethod {
    /// Use the raw values.
    Raw,
    /// log10(x + 1) — FPKM/TPM.
    Log10p1,
    /// ln(x) — metabolomics / proteomics abundance.
    Log,
    /// log2(x + 1) — raw count (RNA-seq), basic transform.
    Log2p1,
    /// log10(x) — peak area (metabolomics).
    Log10,
    /// log2(x) — protein abundance.
    Log2,
}

/// R's `median` (average of the two middle values for even length).
pub fn median(x: &[f64]) -> f64 {
    let mut v = x.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 0 {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    } else {
        v[n / 2]
    }
}

/// R's `mad(x)` (median absolute deviation with constant 1.4826).
pub fn mad(x: &[f64]) -> f64 {
    let m = median(x);
    let mut abs_dev: Vec<f64> = x.iter().map(|&v| (v - m).abs()).collect();
    1.4826 * median(&mut abs_dev)
}

/// R's `var(x)` (sample variance, n-1 denominator).
pub fn variance(x: &[f64]) -> f64 {
    let n = x.len();
    if n < 2 {
        return f64::NAN;
    }
    let m = x.iter().sum::<f64>() / n as f64;
    x.iter().map(|&v| (v - m) * (v - m)).sum::<f64>() / (n as f64 - 1.0)
}

/// Apply a normalization transform.
pub fn normalize(expr: &Mat<f64>, method: NormalizeMethod) -> Mat<f64> {
    let (g, n) = (expr.nrows(), expr.ncols());
    let mut out = Mat::<f64>::zeros(g, n);
    for i in 0..g {
        for j in 0..n {
            let v = expr[(i, j)];
            out[(i, j)] = match method {
                NormalizeMethod::Raw => v,
                NormalizeMethod::Log10p1 => (v + 1.0).log10(),
                NormalizeMethod::Log => v.ln(),
                NormalizeMethod::Log2p1 => (v + 1.0).log2(),
                NormalizeMethod::Log10 => v.log10(),
                NormalizeMethod::Log2 => v.log2(),
            };
        }
    }
    out
}

/// Low-expression filter: keep genes whose expression is below `cutoff` in at
/// most `sample_percentage` fraction of samples (WGCNA FAQ recommendation).
///
/// Returns a keep-mask (one bool per gene).
pub fn low_expression_mask(expr: &Mat<f64>, cutoff: f64, sample_percentage: f64) -> Vec<bool> {
    let n = expr.ncols();
    let max_below = sample_percentage * n as f64;
    (0..expr.nrows())
        .map(|i| {
            let cnt = (0..n).filter(|&j| expr[(i, j)] < cutoff).count();
            (cnt as f64) <= max_below
        })
        .collect()
}

/// Gene-selection method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterMethod {
    Mad,
    Var,
}

/// Select the top `n_keep` genes by MAD or variance. Returns their row indices,
/// sorted by descending score (ties broken by lower index for determinism).
pub fn top_genes_by(expr: &Mat<f64>, n_keep: usize, method: FilterMethod) -> Vec<usize> {
    let g = expr.nrows();
    let n = expr.ncols();
    let mut scores: Vec<(f64, usize)> = (0..g)
        .map(|i| {
            let row: Vec<f64> = (0..n).map(|j| expr[(i, j)]).collect();
            let s = match method {
                FilterMethod::Mad => mad(&row),
                FilterMethod::Var => variance(&row),
            };
            (s, i)
        })
        // 排除零方差（常数）基因：其相关性为 NaN，会导致下游 panic。
        .filter(|(s, _)| s.is_finite() && *s > 0.0)
        .collect();
    // sort by score descending; tie -> lower index first
    scores.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)));
    scores.into_iter().take(n_keep).map(|(_, i)| i).collect()
}

/// Apply a keep-mask to produce a filtered matrix (plus the surviving gene
/// names if provided).
pub fn apply_mask(
    expr: &Mat<f64>,
    keep: &[bool],
    genes: &[String],
) -> Result<(Mat<f64>, Vec<String>)> {
    if keep.len() != expr.nrows() {
        return Err(WgcnaError::InvalidInput(
            "keep mask length != number of genes".into(),
        ));
    }
    let idx: Vec<usize> = (0..expr.nrows()).filter(|&i| keep[i]).collect();
    let g = idx.len();
    let n = expr.ncols();
    let mut out = Mat::<f64>::zeros(g, n);
    for (r, &i) in idx.iter().enumerate() {
        for j in 0..n {
            out[(r, j)] = expr[(i, j)];
        }
    }
    let out_genes: Vec<String> = idx.iter().map(|&i| genes[i].clone()).collect();
    Ok((out, out_genes))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_even_odd() {
        assert!((median(&[1.0, 2.0, 3.0]) - 2.0).abs() < 1e-15);
        assert!((median(&[1.0, 2.0, 3.0, 4.0]) - 2.5).abs() < 1e-15);
    }

    #[test]
    fn variance_known() {
        // var(c(1,2,3,4,5)) = 2.5
        assert!((variance(&[1.0, 2.0, 3.0, 4.0, 5.0]) - 2.5).abs() < 1e-12);
    }

    #[test]
    fn low_expression_mask_keeps_and_removes() {
        let m = faer::mat![[10.0, 10.0, 10.0], [1.0, 1.0, 1.0], [10.0, 1.0, 10.0]];
        // cutoff 5, sample_percentage 0.5 (max_below = 1.5):
        // gene0: 0 below -> keep; gene1: 3 below -> remove; gene2: 1 below -> keep
        let keep = low_expression_mask(&m, 5.0, 0.5);
        assert_eq!(keep, vec![true, false, true]);
    }
}
