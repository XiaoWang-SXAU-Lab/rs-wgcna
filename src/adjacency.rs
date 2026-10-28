//! Adjacency transform, matching WGCNA's `adjacency()`.

use faer::Mat;
use rayon::prelude::*;

use crate::network::{NetworkType, Power};

/// Apply the network-type adjacency transform to a correlation matrix.
///
/// Faithful port of the four branches in WGCNA's C `adjacency()`:
/// - unsigned: `|cor|^power`
/// - signed: `((1+cor)/2)^power`
/// - signed hybrid: `cor^power` if `cor > 0` else `0`
///
/// The diagonal stays `1` (cor(i,i)=1), exactly as in WGCNA.
///
/// Uses `powi` for exactly integral powers and `powf` for real powers.
/// The element-wise transform is parallelized.
pub fn adjacency(cor: &Mat<f64>, power: f64, network_type: NetworkType) -> Mat<f64> {
    let g = cor.nrows();
    let p = Power::new(power);
    let n = g * g;
    // 列优先的平铺索引：idx = i + j*g（faer 默认列主序）
    let data: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|idx| {
            let i = idx % g;
            let j = idx / g;
            let c = cor[(i, j)];
            match network_type {
                NetworkType::Unsigned => p.apply(c.abs()),
                NetworkType::Signed => p.apply((1.0 + c) / 2.0),
                NetworkType::SignedHybrid => {
                    if c > 0.0 {
                        p.apply(c)
                    } else {
                        0.0
                    }
                }
            }
        })
        .collect();
    let mut a = Mat::<f64>::zeros(g, g);
    for j in 0..g {
        for i in 0..g {
            a[(i, j)] = data[i + j * g];
        }
    }
    a
}
