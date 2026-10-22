//! Soft threshold selection — faithful port of WGCNA's `scaleFreeFitIndex`
//! and `pickSoftThreshold` (R/networkFunctions.R, R/Functions.R).

use faer::Mat;

use crate::adjacency::adjacency;
use crate::network::NetworkType;

/// Scale-free topology fit index for a connectivity vector `k`.
///
/// Returns `(Rsquared.SFT, slope.SFT, truncatedExponentialAdjRsquared)`,
/// replicating `scaleFreeFitIndex(k, nBreaks, removeFirst)`.
pub fn scale_free_fit_index(k: &[f64], n_breaks: usize, remove_first: bool) -> (f64, f64, f64) {
    let n = k.len();
    let mut min_k = f64::INFINITY;
    let mut max_k = f64::NEG_INFINITY;
    for &v in k {
        if v < min_k {
            min_k = v;
        }
        if v > max_k {
            max_k = v;
        }
    }
    let dx = max_k - min_k;

    // R 4.5 `cut(k, n)`: interior breaks = seq(min, max, length=n+1),
    // then the two endpoints are extended by 0.1% of the range.
    let cut_lo = min_k - dx / 1000.0;
    let cut_hi = max_k + dx / 1000.0;
    let step = dx / (n_breaks as f64);
    let mut cut_breaks = vec![0.0f64; n_breaks + 1];
    cut_breaks[0] = cut_lo;
    for i in 1..n_breaks {
        cut_breaks[i] = min_k + (i as f64) * step;
    }
    cut_breaks[n_breaks] = cut_hi;

    // Assign each value to its interval (right=TRUE, include.lowest=FALSE):
    // interval i = (cut_breaks[i], cut_breaks[i+1]].
    let mut sum = vec![0.0f64; n_breaks];
    let mut cnt = vec![0usize; n_breaks];
    for &v in k {
        let mut idx = cut_breaks.partition_point(|&b| b < v) as isize - 1;
        if idx < 0 {
            idx = 0;
        }
        if idx >= n_breaks as isize {
            idx = n_breaks as isize - 1;
        }
        let i = idx as usize;
        sum[i] += v;
        cnt[i] += 1;
    }

    let mut dk = vec![0.0f64; n_breaks];
    let mut p_dk = vec![0.0f64; n_breaks];
    for i in 0..n_breaks {
        dk[i] = if cnt[i] == 0 {
            f64::NAN
        } else {
            sum[i] / cnt[i] as f64
        };
        p_dk[i] = cnt[i] as f64 / n as f64;
    }

    // `hist(k, breaks = seq(min,max,length=nBreaks+1), right=TRUE)`: mids.
    // break_i = min + i * (dx/n); mids_i = min + (i + 0.5) * (dx/n).
    let hist_step = dx / (n_breaks as f64);
    let dk2: Vec<f64> = (0..n_breaks)
        .map(|i| min_k + (i as f64 + 0.5) * hist_step)
        .collect();

    for i in 0..n_breaks {
        if dk[i].is_nan() || dk[i] == 0.0 {
            dk[i] = dk2[i];
        }
        if p_dk[i].is_nan() {
            p_dk[i] = 0.0;
        }
    }

    let mut log_dk: Vec<f64> = dk.iter().map(|&d| d.log10()).collect();
    let mut log_p_dk: Vec<f64> = p_dk.iter().map(|&p| (p + 1e-9).log10()).collect();

    if remove_first {
        log_dk.remove(0);
        log_p_dk.remove(0);
    }

    let (slope, r2) = simple_lm(&log_dk, &log_p_dk);
    let z: Vec<f64> = log_dk.iter().map(|&x| 10f64.powf(x)).collect();
    let adj_r2 = multiple_adj_r2(&log_dk, &z, &log_p_dk);

    (r2, slope, adj_r2)
}

/// Simple linear regression y ~ x: returns (slope, R-squared).
fn simple_lm(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len();
    let mx = x.iter().sum::<f64>() / n as f64;
    let my = y.iter().sum::<f64>() / n as f64;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    let mut sxy = 0.0;
    for i in 0..n {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    let slope = if sxx == 0.0 { f64::NAN } else { sxy / sxx };
    let r2 = if sxx == 0.0 || syy == 0.0 {
        f64::NAN
    } else {
        (sxy * sxy) / (sxx * syy)
    };
    (slope, r2)
}

/// Adjusted R-squared of y ~ x1 + x2 (with intercept), via 3x3 normal equations.
fn multiple_adj_r2(x1: &[f64], x2: &[f64], y: &[f64]) -> f64 {
    let n = x1.len();
    let mut a = [[0.0f64; 3]; 3];
    let mut b = [0.0f64; 3];
    for i in 0..n {
        let row = [1.0, x1[i], x2[i]];
        for p in 0..3 {
            for q in 0..3 {
                a[p][q] += row[p] * row[q];
            }
            b[p] += row[p] * y[i];
        }
    }
    let beta = solve3(a, b);
    let my = y.iter().sum::<f64>() / n as f64;
    let mut ss_res = 0.0;
    let mut ss_tot = 0.0;
    for i in 0..n {
        let fit = beta[0] + beta[1] * x1[i] + beta[2] * x2[i];
        ss_res += (y[i] - fit) * (y[i] - fit);
        ss_tot += (y[i] - my) * (y[i] - my);
    }
    let r2 = if ss_tot == 0.0 { f64::NAN } else { 1.0 - ss_res / ss_tot };
    let p = 2.0;
    if r2.is_nan() {
        return r2;
    }
    1.0 - (1.0 - r2) * (n as f64 - 1.0) / (n as f64 - p - 1.0)
}

fn solve3(mut a: [[f64; 3]; 3], mut b: [f64; 3]) -> [f64; 3] {
    for col in 0..3 {
        let mut piv = col;
        for r in (col + 1)..3 {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let d = a[col][col];
        for r in (col + 1)..3 {
            let f = a[r][col] / d;
            for c in col..3 {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = [0.0f64; 3];
    for r in (0..3).rev() {
        let mut s = b[r];
        for c in (r + 1)..3 {
            s -= a[r][c] * x[c];
        }
        x[r] = s / a[r][r];
    }
    x
}

/// Connectivity vector from an adjacency matrix: `k_i = sum_j A_ij - 1`
/// (the adjacency diagonal is 1).
pub fn connectivity(a: &Mat<f64>) -> Vec<f64> {
    let g = a.nrows();
    (0..g)
        .map(|i| {
            let mut s = 0.0;
            for j in 0..g {
                s += a[(i, j)];
            }
            s - 1.0
        })
        .collect()
}

/// Soft threshold selection (port of `pickSoftThreshold`).
///
/// Returns `(power_estimate, rows)` where each row is
/// `[power, R.sq, slope, truncated.R.sq, mean_k, median_k, max_k]`.
pub fn pick_soft_threshold(
    cor: &Mat<f64>,
    network_type: NetworkType,
    power_vector: &[f64],
    rsquared_cut: f64,
    n_breaks: usize,
    remove_first: bool,
) -> (Option<f64>, Vec<[f64; 7]>) {
    let mut rows = Vec::with_capacity(power_vector.len());
    let mut power_estimate = None;
    for &power in power_vector {
        let a = adjacency(cor, power, network_type);
        let k = connectivity(&a);
        let (r2, slope, trunc_r2) = scale_free_fit_index(&k, n_breaks, remove_first);
        let mean_k = k.iter().sum::<f64>() / k.len() as f64;
        let mut ks: Vec<f64> = k.iter().copied().filter(|v| v.is_finite()).collect();
        if ks.is_empty() {
            ks = k.clone();
        }
        ks.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median_k = if ks.len() % 2 == 0 {
            (ks[ks.len() / 2 - 1] + ks[ks.len() / 2]) / 2.0
        } else {
            ks[ks.len() / 2]
        };
        let max_k = ks[ks.len() - 1];
        if power_estimate.is_none() && r2 > rsquared_cut {
            power_estimate = Some(power);
        }
        rows.push([power, r2, slope, trunc_r2, mean_k, median_k, max_k]);
    }
    (power_estimate, rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn connectivity_diag_one() {
        // 2x2 adjacency with diagonal 1: k = 1 - 1 = 0 each.
        let a = faer::mat![[1.0, 0.5], [0.5, 1.0]];
        let k = connectivity(&a);
        assert!((k[0] - 0.5).abs() < 1e-15);
        assert!((k[1] - 0.5).abs() < 1e-15);
    }

    #[test]
    fn simple_lm_slope_r2() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [2.0, 4.0, 6.0, 8.0, 10.0];
        let (slope, r2) = simple_lm(&x, &y);
        assert!((slope - 2.0).abs() < 1e-12);
        assert!((r2 - 1.0).abs() < 1e-12);
    }
}
