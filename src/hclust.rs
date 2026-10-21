//! Hierarchical clustering (average linkage / UPGMA).
//!
//! Faithful port of R's `hclust` (src/library/stats/src/hclust.f, the Murtagh
//! nearest-neighbour algorithm) so that the merge sequence, heights, and leaf
//! order match `stats::hclust(method="average")` / `fastcluster::hclust` exactly.

use faer::Mat;

/// Result of hierarchical clustering, in R's `hclust` convention.
#[derive(Debug, Clone)]
pub struct HClust {
    /// n-1 merges; each entry `(a, b)` with negative = leaf `-(i+1)` and
    /// non-negative = cluster index (1-based).
    pub merge: Vec<(i32, i32)>,
    /// Merge heights, length n-1.
    pub height: Vec<f64>,
    /// Leaf order (1-based leaf indices), length n.
    pub order: Vec<i32>,
}

/// Average-linkage (UPGMA) clustering of a dissimilarity matrix `diss` (n x n,
/// symmetric, diagonal ignored).
pub fn hclust_average(diss: &Mat<f64>) -> HClust {
    let n = diss.nrows();
    if n == 0 {
        return HClust {
            merge: Vec::new(),
            height: Vec::new(),
            order: Vec::new(),
        };
    }
    if n == 1 {
        return HClust {
            merge: Vec::new(),
            height: Vec::new(),
            order: vec![1],
        };
    }

    // Working copy of the dissimilarity matrix (updated in place).
    let mut d = diss.clone();

    // 1-based arrays (index 0 unused) to mirror the Fortran closely.
    let mut flag = vec![true; n + 1];
    let mut membr = vec![0.0f64; n + 1];
    for i in 1..=n {
        membr[i] = 1.0;
    }
    let mut nn = vec![0usize; n + 1];
    let mut disnn = vec![f64::INFINITY; n + 1];

    let mut ia = vec![0i32; n + 1];
    let mut ib = vec![0i32; n + 1];
    let mut crit = vec![0.0f64; n + 1];

    // Initial nearest neighbours (to the right: J > I).
    for i in 1..n {
        let mut dmin = f64::INFINITY;
        let mut jm = 0usize;
        for j in (i + 1)..=n {
            let v = d[(i - 1, j - 1)];
            if dmin > v {
                dmin = v;
                jm = j;
            }
        }
        nn[i] = jm;
        disnn[i] = dmin;
    }

    let mut ncl = n;
    for _ in 0..(n - 1) {
        // Find the active cluster with the smallest nearest-neighbour distance.
        let mut dmin = f64::INFINITY;
        let mut im = 0usize;
        let mut jm = 0usize;
        for i in 1..n {
            if flag[i] && disnn[i] < dmin {
                dmin = disnn[i];
                im = i;
                jm = nn[i];
            }
        }
        ncl -= 1;
        let i2 = im.min(jm);
        let j2 = im.max(jm);
        // Fortran stores at index N-NCL (first merge -> index 1).
        let idx = n - ncl;
        ia[idx] = i2 as i32;
        ib[idx] = j2 as i32;
        crit[idx] = dmin;
        flag[j2] = false;

        // Update dissimilarities from the merged cluster i2 to every other
        // active cluster k (average linkage).
        let mut dmin2 = f64::INFINITY;
        let mut jj = 0usize;
        for k in 1..=n {
            if flag[k] && k != i2 {
                let d1 = d[(i2 - 1, k - 1)];
                let d2 = d[(j2 - 1, k - 1)];
                let newd = (membr[i2] * d1 + membr[j2] * d2) / (membr[i2] + membr[j2]);
                d[(i2 - 1, k - 1)] = newd;
                d[(k - 1, i2 - 1)] = newd;
                if i2 < k {
                    if newd < dmin2 {
                        dmin2 = newd;
                        jj = k;
                    }
                } else if newd < disnn[k] {
                    disnn[k] = newd;
                    nn[k] = i2;
                }
            }
        }
        membr[i2] += membr[j2];
        disnn[i2] = dmin2;
        nn[i2] = jj;

        // Recompute nearest neighbours whose target was i2 or j2.
        for i in 1..n {
            if flag[i] && (nn[i] == i2 || nn[i] == j2) {
                let mut dmin3 = f64::INFINITY;
                let mut jj3 = 0usize;
                for j in (i + 1)..=n {
                    if flag[j] {
                        let v = d[(i - 1, j - 1)];
                        if v < dmin3 {
                            dmin3 = v;
                            jj3 = j;
                        }
                    }
                }
                nn[i] = jj3;
                disnn[i] = dmin3;
            }
        }
    }

    // --- HCASS2: convert IA/IB (representative indices) into R's merge encoding.
    let mut iia = ia.clone();
    let mut iib = ib.clone();
    for i in 1..=(n - 2) {
        let k = ia[i].min(ib[i]);
        for j in (i + 1)..=(n - 1) {
            if ia[j] == k {
                iia[j] = -(i as i32);
            }
            if ib[j] == k {
                iib[j] = -(i as i32);
            }
        }
    }
    for i in 1..=(n - 1) {
        iia[i] = -iia[i];
        iib[i] = -iib[i];
    }
    for i in 1..=(n - 1) {
        if iia[i] > 0 && iib[i] < 0 {
            std::mem::swap(&mut iia[i], &mut iib[i]);
        }
        if iia[i] > 0 && iib[i] > 0 {
            let a = iia[i].min(iib[i]);
            let b = iia[i].max(iib[i]);
            iia[i] = a;
            iib[i] = b;
        }
    }

    let mut merge = Vec::with_capacity(n - 1);
    let mut height = Vec::with_capacity(n - 1);
    for i in 1..=(n - 1) {
        merge.push((iia[i], iib[i]));
        height.push(crit[i]);
    }

    // --- leaf order (IORDER in hclust.f).
    let mut iorder = vec![0i32; n + 1];
    iorder[1] = iia[n - 1];
    iorder[2] = iib[n - 1];
    let mut loc = 2usize;
    for i in (1..=(n - 2)).rev() {
        for j in 1..=loc {
            if iorder[j] == i as i32 {
                iorder[j] = iia[i];
                if j == loc {
                    loc += 1;
                    iorder[loc] = iib[i];
                } else {
                    loc += 1;
                    for k in (j + 2..=loc).rev() {
                        iorder[k] = iorder[k - 1];
                    }
                    iorder[j + 1] = iib[i];
                }
                break;
            }
        }
    }
    let mut order = Vec::with_capacity(n);
    for i in 1..=n {
        order.push(-iorder[i]);
    }

    HClust {
        merge,
        height,
        order,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_points() {
        let d = faer::mat![[0.0, 1.5], [1.5, 0.0]];
        let h = hclust_average(&d);
        assert_eq!(h.merge, vec![(-1, -2)]);
        assert!((h.height[0] - 1.5).abs() < 1e-15);
        assert_eq!(h.order, vec![1, 2]);
    }

    #[test]
    fn three_points_average() {
        // distances: d12=1, d13=3, d23=4
        let d = faer::mat![[0.0, 1.0, 3.0], [1.0, 0.0, 4.0], [3.0, 4.0, 0.0]];
        let h = hclust_average(&d);
        // merge 1&2 (height 1), then {1,2} with 3: avg=(3+4)/2=3.5
        assert_eq!(h.merge[0], (-1, -2));
        assert!((h.height[0] - 1.0).abs() < 1e-15);
        assert_eq!(h.merge[1], (-3, 1));
        assert!((h.height[1] - 3.5).abs() < 1e-15);
        assert_eq!(h.order, vec![3, 1, 2]);
    }
}
