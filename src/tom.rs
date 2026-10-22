//! Topological overlap matrix (TOM).
//!
//! Faithful port of WGCNA's C `tomSimilarity` + `tomSimilarityFromAdj`
//! (src/networkFunctions.c). The entry point takes a *correlation* matrix and
//! reproduces the exact adjacency-type mapping that happens inside
//! `tomSimilarity`, including the `UnsignedKeepSign` adjustment used when
//! `TOMType = "signed"` and `networkType = "unsigned"`.

use faer::linalg::matmul::matmul;
use faer::{Accum, Mat, Par};

use crate::network::{NetworkType, TomDenom, TomType};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AdjType {
    Unsigned,
    UnsignedKeepSign,
    Signed,
    Hybrid,
}

#[inline]
fn adj_transform(c: f64, power: f64, t: AdjType) -> f64 {
    let p = power as i32;
    match t {
        AdjType::Unsigned => c.abs().powi(p),
        AdjType::UnsignedKeepSign => {
            let s = if c < 0.0 { -1.0 } else { 1.0 };
            s * c.abs().powi(p)
        }
        AdjType::Signed => ((1.0 + c) / 2.0).powi(p),
        AdjType::Hybrid => {
            if c > 0.0 {
                c.powi(p)
            } else {
                0.0
            }
        }
    }
}

/// Compute TOM from a correlation matrix.
///
/// `tom_type = TomType::None` returns the adjacency matrix directly (as in WGCNA).
pub fn tom_from_cor(
    cor: &Mat<f64>,
    power: f64,
    network_type: NetworkType,
    tom_type: TomType,
    denom: TomDenom,
    threads: usize,
) -> Mat<f64> {
    let g = cor.nrows();

    // Map network type to adjacency type.
    let mut adj_type = match network_type {
        NetworkType::Unsigned => AdjType::Unsigned,
        NetworkType::Signed => AdjType::Signed,
        NetworkType::SignedHybrid => AdjType::Hybrid,
    };

    // The adjType adjustments from WGCNA tomSimilarity().
    match tom_type {
        TomType::Signed => {
            if adj_type == AdjType::Unsigned {
                adj_type = AdjType::UnsignedKeepSign;
            }
        }
        TomType::SignedNowick | TomType::SignedNowick2 => {
            adj_type = AdjType::UnsignedKeepSign;
        }
        TomType::Unsigned | TomType::Unsigned2 => {
            if adj_type == AdjType::UnsignedKeepSign {
                adj_type = AdjType::Unsigned;
            }
        }
        _ => {}
    }

    // Build adjacency.
    let mut a = Mat::<f64>::zeros(g, g);
    for j in 0..g {
        for i in 0..g {
            a[(i, j)] = adj_transform(cor[(i, j)], power, adj_type);
        }
    }

    // `TomType::None` returns the adjacency directly.
    if tom_type == TomType::None {
        return a;
    }

    // Set diagonal to 1 and compute connectivity (sum of |A| over the row).
    let mut conn = vec![0.0f64; g];
    for i in 0..g {
        a[(i, i)] = 1.0;
        let mut s = 0.0;
        for j in 0..g {
            s += a[(i, j)].abs();
        }
        conn[i] = s;
    }

    // B = A * A (symmetric).
    let mut b = Mat::<f64>::zeros(g, g);
    let par = Par::rayon(threads);
    matmul(&mut b, Accum::Replace, &a, &a, 1.0, par);

    // The "2" variants use the corrected normalization.
    let form = matches!(
        tom_type,
        TomType::Unsigned2 | TomType::Signed2 | TomType::SignedNowick2
    );

    let mut tom = Mat::<f64>::zeros(g, g);
    for j in 0..g {
        for i in (j + 1)..g {
            let aij = a[(i, j)];
            let bij = b[(i, j)];
            let den1 = match denom {
                TomDenom::Min => conn[i].min(conn[j]),
                TomDenom::Mean => (conn[i] + conn[j]) / 2.0,
            };
            let v = match tom_type {
                TomType::Unsigned | TomType::Unsigned2 => {
                    let den = den1 - aij;
                    if den == 0.0 {
                        0.0
                    } else if form {
                        let r = if den <= 1.0 {
                            0.0
                        } else {
                            (bij - 2.0 * aij) / (den - 1.0)
                        };
                        (aij + r) / 2.0
                    } else {
                        (bij - aij) / den
                    }
                }
                TomType::Signed | TomType::Signed2 => {
                    let den = den1 - aij.abs();
                    if den == 0.0 {
                        0.0
                    } else if form {
                        let r = if den <= 1.0 {
                            0.0
                        } else {
                            (bij - 2.0 * aij) / (den - 1.0)
                        };
                        (aij + r).abs() / 2.0
                    } else {
                        (bij - aij).abs() / den
                    }
                }
                TomType::SignedNowick | TomType::SignedNowick2 => {
                    let den = den1 - aij.abs();
                    if den == 0.0 {
                        0.0
                    } else if form {
                        let r = if den <= 1.0 {
                            0.0
                        } else {
                            (bij - 2.0 * aij) / (den - 1.0)
                        };
                        (aij + r) / 2.0
                    } else {
                        (bij - aij) / den
                    }
                }
                TomType::None => aij, // unreachable
            };
            tom[(i, j)] = v;
            tom[(j, i)] = v;
        }
    }

    for i in 0..g {
        tom[(i, i)] = 1.0;
    }

    tom
}
