//! Clustering for the "module overview" comparison: KMeans (Lloyd + k-means++)
//! and fuzzy c-means (Bezdek), plus validity metrics (elbow WSS for kmeans,
//! Xie-Beni index + partition coefficient for c-means). Pure Rust — R is only
//! used for plotting.

use faer::Mat;

/// Seeded RNG (splitmix64) for reproducible initializations.
struct SplitMix64(u64);
impl SplitMix64 {
    fn new(seed: u64) -> Self {
        SplitMix64(seed)
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn rows_of(data: &Mat<f64>) -> Vec<Vec<f64>> {
    (0..data.nrows())
        .map(|i| (0..data.ncols()).map(|j| data[(i, j)]).collect())
        .collect()
}

fn euclid2(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum()
}

/// K-means++ initialization: pick k centers with probability proportional to
/// squared distance to the nearest already-chosen center.
fn kmeans_plusplus(rows: &[Vec<f64>], k: usize, rng: &mut SplitMix64) -> Vec<Vec<f64>> {
    let g = rows.len();
    let mut centers = Vec::with_capacity(k);
    let first = (rng.next_f64() * g as f64) as usize % g;
    centers.push(rows[first].clone());
    while centers.len() < k {
        let d2: Vec<f64> = rows
            .iter()
            .map(|r| {
                centers
                    .iter()
                    .map(|c| euclid2(r, c))
                    .fold(f64::INFINITY, f64::min)
            })
            .collect();
        let total: f64 = d2.iter().sum();
        if total <= 0.0 {
            // all points are centers already (degenerate); pick any remaining
            let next = (rng.next_f64() * g as f64) as usize % g;
            centers.push(rows[next].clone());
            continue;
        }
        let mut t = rng.next_f64() * total;
        let mut pick = g - 1;
        for (i, &d) in d2.iter().enumerate() {
            if t < d {
                pick = i;
                break;
            }
            t -= d;
        }
        centers.push(rows[pick].clone());
    }
    centers
}

/// One Lloyd's run. Returns (assignment, centers, wss).
fn kmeans_once(
    rows: &[Vec<f64>],
    k: usize,
    rng: &mut SplitMix64,
    max_iter: usize,
) -> (Vec<usize>, Vec<Vec<f64>>, f64) {
    let g = rows.len();
    let n = rows[0].len();
    let mut centers = kmeans_plusplus(rows, k, rng);
    let mut assign = vec![0usize; g];
    for _ in 0..max_iter {
        let mut changed = false;
        for i in 0..g {
            let mut best = usize::MAX;
            let mut bd = f64::INFINITY;
            for (c, center) in centers.iter().enumerate() {
                let d = euclid2(&rows[i], center);
                if d < bd {
                    bd = d;
                    best = c;
                }
            }
            if assign[i] != best {
                assign[i] = best;
                changed = true;
            }
        }
        // recompute centers
        let mut sums = vec![vec![0.0; n]; k];
        let mut counts = vec![0usize; k];
        for i in 0..g {
            let c = assign[i];
            counts[c] += 1;
            for j in 0..n {
                sums[c][j] += rows[i][j];
            }
        }
        for c in 0..k {
            if counts[c] > 0 {
                for j in 0..n {
                    centers[c][j] = sums[c][j] / counts[c] as f64;
                }
            }
        }
        if !changed {
            break;
        }
    }
    let wss: f64 = (0..g)
        .map(|i| euclid2(&rows[i], &centers[assign[i]]))
        .sum();
    (assign, centers, wss)
}

/// KMeans with multiple restarts (best WSS). Returns (assignment, centers, wss).
pub fn kmeans(
    data: &Mat<f64>,
    k: usize,
    seed: u64,
    max_iter: usize,
    nstart: usize,
) -> (Vec<usize>, Vec<Vec<f64>>, f64) {
    let rows = rows_of(data);
    let mut rng = SplitMix64::new(seed);
    let mut best: Option<(Vec<usize>, Vec<Vec<f64>>, f64)> = None;
    for _ in 0..nstart {
        let r = kmeans_once(&rows, k, &mut rng, max_iter);
        if best.as_ref().map_or(true, |(_, _, bw)| r.2 < *bw) {
            best = Some(r);
        }
    }
    best.expect("kmeans nstart > 0")
}

/// Fuzzy c-means (Bezdek). Returns (membership [g × k], centers [k × n]).
pub fn cmeans(
    data: &Mat<f64>,
    k: usize,
    seed: u64,
    max_iter: usize,
    fuzziness: f64,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let rows = rows_of(data);
    let g = rows.len();
    let n = rows[0].len();
    let mut rng = SplitMix64::new(seed);
    let m = fuzziness.max(1.1);

    // Center-based init (k-means++) avoids the degenerate uniform-membership
    // fixed point that random membership init can fall into.
    let mut centers = kmeans_plusplus(&rows, k, &mut rng);

    let compute_u = |rows: &[Vec<f64>], centers: &[Vec<f64>]| -> Vec<Vec<f64>> {
        let exp = 1.0 / (m - 1.0);
        rows.iter()
            .map(|r| {
                let d2: Vec<f64> = centers.iter().map(|c| euclid2(r, c)).collect();
                let mut new = vec![0.0; k];
                let zero = d2.iter().position(|&d| d < 1e-12);
                if let Some(zc) = zero {
                    new[zc] = 1.0;
                } else {
                    for c in 0..k {
                        let mut s = 0.0;
                        for cc in 0..k {
                            s += (d2[c] / d2[cc]).powf(exp);
                        }
                        new[c] = 1.0 / s;
                    }
                }
                new
            })
            .collect()
    };

    let mut u = compute_u(&rows, &centers);
    for _ in 0..max_iter {
        // update centers
        for c in 0..k {
            let mut denom = 0.0;
            let mut num = vec![0.0; n];
            for i in 0..g {
                let w = u[i][c].powf(m);
                denom += w;
                for j in 0..n {
                    num[j] += w * rows[i][j];
                }
            }
            for j in 0..n {
                centers[c][j] = if denom > 0.0 { num[j] / denom } else { 0.0 };
            }
        }
        // update membership
        let new_u = compute_u(&rows, &centers);
        let mut max_delta = 0.0f64;
        for i in 0..g {
            for c in 0..k {
                let delta = (new_u[i][c] - u[i][c]).abs();
                if delta > max_delta {
                    max_delta = delta;
                }
            }
        }
        u = new_u;
        if max_delta < 1e-4 {
            break;
        }
    }
    (u, centers)
}

/// Elbow for kmeans: WSS over k = 2..kmax. Returns (wss, ks).
pub fn kmeans_wss(data: &Mat<f64>, kmax: usize, seed: u64) -> (Vec<f64>, Vec<usize>) {
    let ks: Vec<usize> = (2..=kmax).collect();
    let wss: Vec<f64> = ks
        .iter()
        .map(|&k| kmeans(data, k, seed, 100, 5).2)
        .collect();
    (wss, ks)
}

/// Xie-Beni index + partition coefficient for fuzzy c-means over k = 2..kmax.
pub fn cmeans_validity(
    data: &Mat<f64>,
    kmax: usize,
    seed: u64,
    fuzziness: f64,
) -> (Vec<f64>, Vec<f64>, Vec<usize>) {
    let ks: Vec<usize> = (2..=kmax).collect();
    let mut xb = Vec::with_capacity(ks.len());
    let mut pc = Vec::with_capacity(ks.len());
    let rows = rows_of(data);
    let g = rows.len();
    for &k in &ks {
        let (u, centers) = cmeans(data, k, seed, 100, fuzziness);
        // partition coefficient
        let pcv: f64 = u.iter().flatten().map(|v| v * v).sum::<f64>() / g as f64;
        // Xie-Beni: compactness / (g * separation)
        let mut compact = 0.0;
        let mut sep = f64::INFINITY;
        for i in 0..g {
            // compactness to own (argmax) center
            let hard = (0..k)
                .max_by(|&a, &b| u[i][a].partial_cmp(&u[i][b]).unwrap())
                .unwrap();
            compact += euclid2(&rows[i], &centers[hard]);
        }
        for a in 0..k {
            for b in (a + 1)..k {
                let d = euclid2(&centers[a], &centers[b]);
                if d < sep {
                    sep = d;
                }
            }
        }
        let xbv = if sep > 0.0 { compact / (g as f64 * sep) } else { f64::INFINITY };
        xb.push(xbv);
        pc.push(pcv);
    }
    (xb, pc, ks)
}

/// Reorder clusters by size (largest first) and relabel to C1..Ck.
pub fn relabel_by_size(assign: &[usize], k: usize) -> Vec<usize> {
    let mut counts = vec![0usize; k];
    for &a in assign {
        counts[a] += 1;
    }
    let mut order: Vec<usize> = (0..k).collect();
    order.sort_by_key(|&c| std::cmp::Reverse(counts[c]));
    let mut map = vec![0usize; k];
    for (new, &old) in order.iter().enumerate() {
        map[old] = new;
    }
    assign.iter().map(|&a| map[a]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cmeans_finds_two_clusters() {
        // 2 clear groups in 4D
        let mut m = Mat::<f64>::zeros(20, 4);
        for i in 0..10 {
            for j in 0..4 {
                m[(i, j)] = 10.0 + (i as f64 * 0.1) + (j as f64);
            }
        }
        for i in 10..20 {
            for j in 0..4 {
                m[(i, j)] = -10.0 - (i as f64 * 0.1) + (j as f64);
            }
        }
        let (u, _centers) = cmeans(&m, 2, 42, 50, 2.0);
        let pc: f64 = u.iter().flatten().map(|v| v * v).sum::<f64>() / 20.0;
        eprintln!("PC = {}", pc);
        assert!(pc > 0.9, "PC should be near 1 for clean clusters, got {}", pc);
    }
}
