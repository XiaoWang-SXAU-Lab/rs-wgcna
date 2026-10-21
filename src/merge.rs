//! mergeCloseModules: merge modules whose eigengenes are highly correlated.
//!
//! Faithful port of WGCNA's `mergeCloseModules` for the single-set, pearson,
//! `useAbs = FALSE`, `relabel = TRUE` path used by `blockwiseModules`.

use std::collections::BTreeMap;

use faer::Mat;

use crate::correlation::pearson;
use crate::eigengene::module_eigengene;
use crate::hclust::hclust_average;

/// First entries of WGCNA's `standardColors()`.
pub const STANDARD_COLORS: &[&str] = &[
    "turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink",
    "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue",
    "lightcyan", "grey60", "lightgreen", "lightyellow", "royalblue", "darkred",
    "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white",
    "skyblue", "saddlebrown", "steelblue", "paleturquoise", "violet",
    "darkolivegreen", "darkmagenta", "sienna3", "yellowgreen", "skyblue3",
    "plum1", "orangered4", "mediumpurple3", "lightsteelblue1", "lightcyan1",
    "ivory", "floralwhite", "darkorange2", "brown4", "bisque4", "darkslateblue",
    "plum2", "thistle2",
];

fn find(parent: &mut [usize], x: usize) -> usize {
    if parent[x] != x {
        parent[x] = find(parent, parent[x]);
    }
    parent[x]
}

/// Partition leaves into connected components by cutting the dendrogram at
/// height `h` (matches `cutree(tree, h = h)` for monotone methods).
fn cutree_partition(merge: &[(i32, i32)], height: &[f64], h: f64) -> Vec<usize> {
    let n = merge.len() + 1;
    let mut parent: Vec<usize> = (0..n).collect();
    let mut cluster_rep: Vec<usize> = vec![0; merge.len()];

    for (i, &(a, b)) in merge.iter().enumerate() {
        if height[i] > h {
            break; // heights are non-decreasing for average linkage
        }
        let la = if a < 0 {
            (-a - 1) as usize
        } else {
            cluster_rep[a as usize - 1]
        };
        let lb = if b < 0 {
            (-b - 1) as usize
        } else {
            cluster_rep[b as usize - 1]
        };
        let ra = find(&mut parent, la);
        let rb = find(&mut parent, lb);
        if ra != rb {
            parent[rb] = ra;
        }
        cluster_rep[i] = find(&mut parent, la);
    }

    (0..n).map(|x| find(&mut parent, x)).collect()
}

/// Merge close modules and relabel by size (grey stays grey).
///
/// `expr` is genes x samples; `colors` is one color per gene (e.g. "grey" for
/// unassigned). Returns merged + relabeled colors.
pub fn merge_close_modules(
    expr: &Mat<f64>,
    colors: &[String],
    cut_height: f64,
    threads: usize,
) -> Vec<String> {
    let n_genes = colors.len();
    let grey = "grey".to_string();
    let mut merged = colors.to_vec();

    loop {
        // Sorted unique non-grey modules.
        let mut mods: Vec<String> = merged.iter().cloned().collect();
        mods.sort();
        mods.dedup();
        let mods: Vec<String> = mods.into_iter().filter(|c| *c != grey).collect();
        if mods.len() < 2 {
            break;
        }

        // Module eigengenes (samples x modules).
        let n_samples = expr.ncols();
        let k = mods.len();
        let mut mes = Mat::<f64>::zeros(n_samples, k);
        for (mi, m) in mods.iter().enumerate() {
            let idx: Vec<usize> = (0..n_genes).filter(|&i| merged[i] == *m).collect();
            let me = module_eigengene(expr, &idx);
            for j in 0..n_samples {
                mes[(j, mi)] = me[j];
            }
        }

        // 1 - cor(MEs), then hierarchical clustering + cut.
        let mes_t = mes.transpose().to_owned();
        let cor_me = pearson(&mes_t, threads).expect("cor of MEs");
        let mut diss = Mat::<f64>::zeros(k, k);
        for i in 0..k {
            for j in 0..k {
                diss[(i, j)] = 1.0 - cor_me[(i, j)];
            }
        }
        for i in 0..k {
            diss[(i, i)] = 0.0;
        }
        let tree = hclust_average(&diss);
        let part = cutree_partition(&tree.merge, &tree.height, cut_height);

        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (mi, &c) in part.iter().enumerate() {
            groups.entry(c).or_default().push(mi);
        }

        let mut did_merge = false;
        for members in groups.values() {
            if members.len() > 1 {
                let first = mods[members[0]].clone();
                for &mi in &members[1..] {
                    let c = mods[mi].clone();
                    for g in 0..n_genes {
                        if merged[g] == c {
                            merged[g] = first.clone();
                        }
                    }
                }
                did_merge = true;
            }
        }

        if !did_merge {
            break;
        }
    }

    relabel(&merged, &grey)
}

fn relabel(colors: &[String], grey: &str) -> Vec<String> {
    let mut counts: BTreeMap<&String, usize> = BTreeMap::new();
    for c in colors {
        *counts.entry(c).or_insert(0) += 1;
    }
    let mut sorted: Vec<&String> = counts
        .keys()
        .copied()
        .filter(|c| c.as_str() != grey)
        .collect();
    sorted.sort_by(|a, b| counts[b].cmp(&counts[a]).then(a.cmp(b)));

    let mut map: BTreeMap<&String, String> = BTreeMap::new();
    for (rank, c) in sorted.iter().enumerate() {
        if rank < STANDARD_COLORS.len() {
            map.insert(c, STANDARD_COLORS[rank].to_string());
        }
    }

    colors
        .iter()
        .map(|c| {
            if c.as_str() == grey {
                grey.to_string()
            } else {
                map.get(c).cloned().unwrap_or_else(|| c.clone())
            }
        })
        .collect()
}
