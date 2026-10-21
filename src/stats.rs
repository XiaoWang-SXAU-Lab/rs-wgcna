//! Statistical helpers: Pearson correlation of two vectors and Student's t CDF
//! (needed for correlation p-values, matching `cor.test`/`corPvalueStudent`).

/// Pearson correlation of two equal-length vectors. Returns NaN if either has
/// zero variance.
pub fn pearson_vec(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len();
    if n < 2 {
        return f64::NAN;
    }
    let ma = a.iter().sum::<f64>() / n as f64;
    let mb = b.iter().sum::<f64>() / n as f64;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    let mut sxy = 0.0;
    for i in 0..n {
        let da = a[i] - ma;
        let db = b[i] - mb;
        sxx += da * da;
        syy += db * db;
        sxy += da * db;
    }
    let den = (sxx * syy).sqrt();
    if den == 0.0 {
        f64::NAN
    } else {
        sxy / den
    }
}

/// Regularized incomplete beta function I_x(a, b) (continued fraction).
fn betai(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let bt = (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b)
        + a * x.ln()
        + b * (1.0 - x).ln())
    .exp();
    // Continued fraction (Numerical Recipes betacf).
    let fpmin = 1.0e-300;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < fpmin {
        d = fpmin;
    }
    d = 1.0 / d;
    let mut h = d;
    let mut m2: f64;
    for m in 1..200 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let mut aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = 1.0 + aa / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = 1.0 + aa / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1.0e-15 {
            break;
        }
    }
    let _ = m2;
    bt * h / a
}

/// Log-gamma via Lanczos approximation.
fn ln_gamma(x: f64) -> f64 {
    if x < 0.5 {
        // reflection formula
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln()
            - ln_gamma(1.0 - x);
    }
    let cof = [
        76.180_091_729_471_46,
        -86.505_320_329_416_77,
        24.014_098_240_830_91,
        -1.231_739_572_450_155,
        0.120_865_097_386_617_9e-2,
        -0.539_523_938_495_3e-5,
    ];
    let x = x - 1.0;
    let tmp = x + 5.5;
    let tmp = tmp - (x + 0.5) * tmp.ln();
    let mut ser = 1.000_000_000_190_015;
    for (i, &c) in cof.iter().enumerate() {
        ser += c / (x + (i as f64) + 1.0);
    }
    -tmp + (2.506_628_274_631_000_5 * ser).ln()
}

/// Two-sided p-value of a Pearson correlation with `n` samples, matching
/// R's `cor.test` / WGCNA `corPvalueStudent`.
pub fn cor_pvalue(r: f64, n: usize) -> f64 {
    if n <= 2 || !r.is_finite() || r.abs() >= 1.0 {
        // Degenerate cases; match R behavior closely enough for parity tests.
        if r.abs() >= 1.0 {
            return 0.0;
        }
        return f64::NAN;
    }
    let df = (n - 2) as f64;
    let t = r * (df / (1.0 - r * r)).sqrt();
    // two-sided: 2 * P(T <= -|t|)
    2.0 * student_t_cdf(-t.abs(), df)
}

/// Student's t cumulative distribution function.
fn student_t_cdf(t: f64, df: f64) -> f64 {
    let x = df / (df + t * t);
    let ib = betai(df / 2.0, 0.5, x);
    if t < 0.0 {
        0.5 * ib
    } else {
        1.0 - 0.5 * ib
    }
}

/// Series expansion of the regularized lower incomplete gamma P(a, x)
/// (Numerical Recipes `gser`), accurate to ~1e-15.
fn gser(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..200 {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * 1e-15 {
            break;
        }
    }
    sum * (-x + a * x.ln() - gln).exp()
}

/// Continued fraction for the regularized upper incomplete gamma Q(a, x)
/// (Numerical Recipes `gcf`), accurate to ~1e-15.
fn gcf(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let fpmin = 1e-300;
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / fpmin;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..200 {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = b + an / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-15 {
            break;
        }
    }
    (-x + a * x.ln() - gln).exp() * h
}

/// Regularized upper incomplete gamma Q(a, x).
fn gammq(a: f64, x: f64) -> f64 {
    if x < a + 1.0 {
        1.0 - gser(a, x)
    } else {
        gcf(a, x)
    }
}

/// Complementary error function, accurate to ~1e-15 (via incomplete gamma).
fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        2.0 - erfc(-x)
    } else {
        gammq(0.5, x * x)
    }
}

/// Standard normal CDF (R's `pnorm`), full double precision.
pub fn normal_cdf(x: f64) -> f64 {
    let y = x / std::f64::consts::SQRT_2;
    if y >= 0.0 {
        1.0 - 0.5 * erfc(y)
    } else {
        0.5 * erfc(-y)
    }
}

/// Fisher z-transformation p-value of a correlation (R's / WGCNA's
/// `corPvalueFisher`).
pub fn cor_pvalue_fisher(cor: f64, n_samples: usize) -> f64 {
    let z = (0.5 * ((1.0 + cor) / (1.0 - cor)).ln() * (n_samples as f64 - 3.0).sqrt()).abs();
    2.0 * normal_cdf(-z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pearson_simple() {
        let a = [1.0, 2.0, 3.0, 4.0, 5.0];
        let b = [2.0, 4.0, 6.0, 8.0, 10.0];
        let r = pearson_vec(&a, &b);
        assert!((r - 1.0).abs() < 1e-12);
    }

    #[test]
    fn pearson_negative() {
        let a = [1.0, 2.0, 3.0, 4.0, 5.0];
        let b = [5.0, 4.0, 3.0, 2.0, 1.0];
        let r = pearson_vec(&a, &b);
        assert!((r + 1.0).abs() < 1e-12);
    }

    #[test]
    fn pvalue_monotonic() {
        // |r| larger -> smaller p.
        assert!(cor_pvalue(0.9, 20) < cor_pvalue(0.5, 20));
        assert!(cor_pvalue(0.5, 20) < cor_pvalue(0.1, 20));
    }

    #[test]
    fn pvalue_perfect() {
        assert_eq!(cor_pvalue(1.0, 10), 0.0);
    }
}
