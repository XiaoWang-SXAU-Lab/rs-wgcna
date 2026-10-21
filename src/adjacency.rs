//! Adjacency transform, matching WGCNA's `adjacency()`.

use faer::Mat;

use crate::network::NetworkType;

/// Apply the network-type adjacency transform to a correlation matrix.
///
/// Faithful port of the four branches in WGCNA's C `adjacency()`:
/// - unsigned: `|cor|^power`
/// - signed: `((1+cor)/2)^power`
/// - signed hybrid: `cor^power` if `cor > 0` else `0`
///
/// The diagonal stays `1` (cor(i,i)=1), exactly as in WGCNA.
pub fn adjacency(cor: &Mat<f64>, power: f64, network_type: NetworkType) -> Mat<f64> {
    let g = cor.nrows();
    let mut a = Mat::<f64>::zeros(g, g);
    for j in 0..g {
        for i in 0..g {
            let c = cor[(i, j)];
            a[(i, j)] = match network_type {
                NetworkType::Unsigned => c.abs().powf(power),
                NetworkType::Signed => ((1.0 + c) / 2.0).powf(power),
                NetworkType::SignedHybrid => {
                    if c > 0.0 {
                        c.powf(power)
                    } else {
                        0.0
                    }
                }
            };
        }
    }
    a
}
