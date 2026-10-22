//! DESeq2 variance-stabilizing transformation (vst), pure Rust — no R needed.
//!
//! Faithful port of DESeq2's `varianceStabilizingTransformation` (blind=TRUE,
//! fitType="parametric"): median-ratio size factors → gene-wise NB dispersion
//! MLE → parametric Gamma-GLM fit (dispersion = a0 + a1/mean) → vst formula.

use faer::Mat;

use crate::stats::ln_gamma;

/// Median ratio size factors (DESeq2 default, type="ratio").
fn size_factors(counts: &Mat<f64>) -> Vec<f64> {
    let g = counts.nrows();
    let n = counts.ncols();
    let mut loggm = vec![0.0f64; g];
    for i in 0..g {
        let mut s = 0.0;
        for j in 0..n {
            s += counts[(i, j)].ln();
        }
        loggm[i] = s / n as f64;
    }
    let mut sf = vec![0.0f64; n];
    for j in 0..n {
        let mut ratios: Vec<f64> = Vec::new();
        for i in 0..g {
            if loggm[i].is_finite() && counts[(i, j)] > 0.0 {
                ratios.push(counts[(i, j)].ln() - loggm[i]);
            }
        }
        ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let m = ratios.len();
        sf[j] = if m == 0 {
            0.0
        } else if m % 2 == 0 {
            ((ratios[m / 2 - 1] + ratios[m / 2]) / 2.0).exp()
        } else {
            ratios[m / 2].exp()
        };
    }
    sf
}

/// Negative NB log-likelihood + Cox-Reid correction (DESeq2 useCR=TRUE).
/// CR term: -0.5 * log(det(X'WX)), intercept-only -> -0.5 * log(sum mu/(1+alpha*mu)).
fn nb_neg_ll(log_alpha: f64, y: &[f64], mu: &[f64]) -> f64 {
    let alpha = log_alpha.exp();
    let r = 1.0 / alpha;
    let mut ll = 0.0;
    let mut cr = 0.0;
    for j in 0..y.len() {
        let yj = y[j];
        let mj = mu[j];
        if mj <= 0.0 {
            continue;
        }
        ll += ln_gamma(yj + r) - ln_gamma(r) + yj * (alpha * mj).ln()
            - (yj + r) * (1.0 + alpha * mj).ln();
        cr += mj / (1.0 + alpha * mj);
    }
    -ll + 0.5 * cr.ln()
}

fn brent_minimize<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, tol: f64) -> f64 {
    let gr = (5.0f64.sqrt() - 1.0) / 2.0;
    let mut lo = a;
    let mut hi = b;
    let mut x = a + gr * (b - a);
    let mut w = x;
    let mut v = x;
    let mut fx = f(x);
    let mut fw = fx;
    let mut fv = fx;
    for _ in 0..100 {
        let m = 0.5 * (lo + hi);
        if (hi - lo).abs() < tol {
            break;
        }
        let r = (w - x) * (fw - fv);
        let q = (v - x) * (fv - fx);
        let p = (v - x) * q - (w - x) * r;
        let q2 = 2.0 * (q - r);
        let mut u = x;
        let mut use_parabola = false;
        if q2.abs() > 1e-12 {
            u = x - p / q2;
            use_parabola = u > lo + tol && u < hi - tol && (u - m).abs() < 0.5 * (hi - lo);
        }
        if !use_parabola {
            u = if x < m { x + gr * (hi - x) } else { x - gr * (x - lo) };
        }
        let fu = f(u);
        if fu <= fx {
            if u >= x {
                lo = x;
            } else {
                hi = x;
            }
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        } else {
            if u < x {
                lo = u;
            } else {
                hi = u;
            }
            if fu <= fw || w == x {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    x
}

fn fit_dispersion(y: &[f64], mu: &[f64]) -> f64 {
    let min_disp: f64 = 1e-8;
    let max_disp: f64 = 10.0f64.max(y.len() as f64);
    let log_alpha = brent_minimize(|t| nb_neg_ll(t, y, mu), min_disp.ln(), max_disp.ln(), 1e-6);
    log_alpha.exp()
}

fn gamma_glm_identity(x: &[f64], y: &[f64], start: (f64, f64)) -> (f64, f64) {
    let n = x.len();
    let (mut a0, mut a1) = start;
    for _ in 0..50 {
        let mut s00 = 0.0;
        let mut s10 = 0.0;
        let mut s11 = 0.0;
        let mut b0 = 0.0;
        let mut b1 = 0.0;
        for i in 0..n {
            let mu = a0 + a1 * x[i];
            if mu <= 0.0 {
                continue;
            }
            let w = 1.0 / (mu * mu);
            s00 += w;
            s10 += w * x[i];
            s11 += w * x[i] * x[i];
            b0 += w * y[i];
            b1 += w * x[i] * y[i];
        }
        let det = s00 * s11 - s10 * s10;
        if det.abs() < 1e-12 {
            break;
        }
        let na0 = (s11 * b0 - s10 * b1) / det;
        let na1 = (s00 * b1 - s10 * b0) / det;
        let change = ((na0 - a0).powi(2) + (na1 - a1).powi(2)).sqrt();
        a0 = na0;
        a1 = na1;
        if change < 1e-9 {
            break;
        }
    }
    (a0, a1)
}

fn parametric_fit(means: &[f64], disps: &[f64]) -> (f64, f64) {
    let (mut a0, mut a1) = (0.1, 1.0);
    for _ in 0..11 {
        let mut x: Vec<f64> = Vec::new();
        let mut y: Vec<f64> = Vec::new();
        for i in 0..means.len() {
            let fitted = a0 + a1 / means[i];
            if fitted <= 0.0 {
                continue;
            }
            let res = disps[i] / fitted;
            if res > 1e-4 && res < 15.0 {
                x.push(1.0 / means[i]);
                y.push(disps[i]);
            }
        }
        if x.is_empty() {
            break;
        }
        let old = (a0, a1);
        (a0, a1) = gamma_glm_identity(&x, &y, (a0, a1));
        if !(a0 > 0.0 && a1 > 0.0) {
            break;
        }
        let rel = (a0 / old.0).ln().powi(2) + (a1 / old.1).ln().powi(2);
        if rel < 1e-6 {
            break;
        }
    }
    (a0, a1)
}

fn vst_value(q: f64, a0: f64, a1: f64) -> f64 {
    let inner = a0 * q * (1.0 + a1 + a0 * q);
    let num = 1.0 + a1 + 2.0 * a0 * q + 2.0 * inner.sqrt();
    (num / (4.0 * a0)).ln() / std::f64::consts::LN_2
}

/// Variance-stabilizing transform of a raw count matrix (genes × samples).
/// Returns the vst-transformed matrix.
pub fn vst_transform(counts: &Mat<f64>) -> Mat<f64> {
    let g = counts.nrows();
    let n = counts.ncols();
    let sf = size_factors(counts);

    let mut q = Mat::<f64>::zeros(g, n);
    for i in 0..g {
        for j in 0..n {
            q[(i, j)] = counts[(i, j)] / sf[j];
        }
    }

    let mut means = vec![0.0f64; g];
    for i in 0..g {
        means[i] = (0..n).map(|j| q[(i, j)]).sum::<f64>() / n as f64;
    }

    let mut disps = vec![0.0f64; g];
    for i in 0..g {
        let y: Vec<f64> = (0..n).map(|j| counts[(i, j)]).collect();
        // intercept-only fitted mean: mu_ij = sf_j * rowMeans(counts/sf)
        let mu: Vec<f64> = sf.iter().map(|s| s * means[i]).collect();
        disps[i] = fit_dispersion(&y, &mu);
    }

    let (a0, a1) = parametric_fit(&means, &disps);

    let mut out = Mat::<f64>::zeros(g, n);
    for i in 0..g {
        for j in 0..n {
            out[(i, j)] = vst_value(q[(i, j)], a0, a1);
        }
    }
    out
}
