//! Dynamic tree cut (the "hybrid" method of `dynamicTreeCut::cutreeDynamic`).
//!
//! Faithful port of `cutreeHybrid` (dynamicTreeCut 1.63) for the parameter
//! combination used by WGCNA `blockwiseModules` defaults: no external branch
//! split functions, `useMedoids = FALSE`, `pamStage = TRUE`,
//! `pamRespectsDendro = TRUE`, `respectSmallClusters = TRUE`.

use faer::Mat;

/// Core size used by the dynamic tree cut (`dynamicTreeCut:::.CoreSize`).
fn core_size(branch_size: usize, min_cluster_size: usize) -> usize {
    let base = min_cluster_size as f64 / 2.0 + 1.0;
    if base < branch_size as f64 {
        (base + (branch_size as f64 - base).sqrt()) as usize
    } else {
        branch_size
    }
}

/// Interpolate a value from `data` at (possibly fractional) `index`
/// (`dynamicTreeCut:::.interpolate`). Uses round-half-to-even like R's `round`.
fn interpolate(data: &[f64], index: f64) -> f64 {
    let i = r_round(index);
    let n = data.len() as f64;
    if i < 1.0 {
        return data[0];
    }
    if i >= n {
        return data[data.len() - 1];
    }
    let r = index - i;
    data[i as usize - 1] * (1.0 - r) + data[i as usize] * r
}

/// R's `round()`: round half to even.
fn r_round(x: f64) -> f64 {
    let f = x.floor();
    let d = x - f;
    if d < 0.5 {
        f
    } else if d > 0.5 {
        f + 1.0
    } else {
        // exactly half: round to even
        if f % 2.0 == 0.0 {
            f
        } else {
            f + 1.0
        }
    }
}

/// R's `rank(x, ties.method = "first")`.
fn rank_first(x: &[f64]) -> Vec<usize> {
    let n = x.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| x[a].partial_cmp(&x[b]).unwrap().then(a.cmp(&b)));
    let mut ranks = vec![0usize; n];
    for (r, &i) in idx.iter().enumerate() {
        ranks[i] = r + 1;
    }
    ranks
}

/// Dynamic tree cut (hybrid). `merge`/`height` are the R-style dendrogram
/// (from `hclust_average`); `dist` is the n x n dissimilarity matrix.
///
/// Returns the cluster labels (0 = unassigned/grey, 1..k = modules in
/// decreasing size order), matching `cutreeDynamic(..., method="hybrid")`.
#[allow(clippy::too_many_arguments)]
pub fn cutree_hybrid(
    merge: &[(i32, i32)],
    height: &[f64],
    dist: &Mat<f64>,
    min_cluster_size: usize,
    deep_split: usize, // cutreeDynamic deepSplit (0..4); +1 applied internally
    cut_height: f64,   // detectCutHeight
    pam_stage: bool,
    pam_respects_dendro: bool,
) -> Vec<i32> {
    let n_merge = height.len();
    let n_points = n_merge + 1;
    if n_merge < 1 {
        return vec![0; n_points];
    }

    // --- setup ---
    let ref_quantile = 0.05;
    let mut ref_merge = (n_merge as f64 * ref_quantile).round() as usize;
    if ref_merge < 1 {
        ref_merge = 1;
    }
    let ref_height = height[ref_merge - 1];

    let mut cut_height = cut_height;
    if cut_height > height.iter().cloned().fold(f64::NEG_INFINITY, f64::max) {
        cut_height = height.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    }
    let max_pam_dist = cut_height;

    let n_merge_below_cut = height.iter().filter(|&&h| h <= cut_height).count();
    if n_merge_below_cut < min_cluster_size {
        return vec![0; n_points];
    }

    let def_mcs = [0.64, 0.73, 0.82, 0.91, 0.95];
    let def_mg: Vec<f64> = def_mcs.iter().map(|&m| (1.0 - m) * 3.0 / 4.0).collect();
    let deep_split = deep_split + 1; // R: deepSplit = deepSplit + 1
    let max_core_scatter = interpolate(&def_mcs, deep_split as f64);
    let min_gap = interpolate(&def_mg, deep_split as f64);
    let max_abs_core_scatter = ref_height + max_core_scatter * (cut_height - ref_height);
    let min_abs_gap = min_gap * (cut_height - ref_height);
    let min_split_height = 0.0;
    let min_abs_split_height = ref_height + min_split_height * (cut_height - ref_height);

    // --- branch state (1-based; index 0 unused) ---
    let mx_branches = n_merge_below_cut;
    let mut b = BranchState::new(mx_branches);

    // IndMergeToBranch (1-based over merges), onBranch (1-based over points).
    let mut ind_merge_to_branch = vec![0usize; n_merge + 1];
    let mut on_branch = vec![0usize; n_points + 1];
    let mut n_branches = 0usize;
    let mut root_branch = 0usize;

    // --- main loop over merges ---
    for m in 1..=n_merge {
        if height[m - 1] > cut_height {
            continue;
        }
        let c1 = merge[m - 1].0;
        let c2 = merge[m - 1].1;

        if c1 < 0 && c2 < 0 {
            // two leaves -> new basic branch
            n_branches += 1;
            b.set_basic(n_branches, (-c1) as usize, (-c2) as usize, height[m - 1]);
            ind_merge_to_branch[m] = n_branches;
            root_branch = n_branches;
        } else if c1.signum() * c2.signum() < 0 {
            // one leaf + one cluster
            let clust = ind_merge_to_branch[c1.max(c2) as usize];
            if clust == 0 {
                panic!("internal error: merge with no associated cluster");
            }
            let gene = (-c1.min(c2)) as usize;
            let ns = b.n_singletons[clust] + 1;
            let nm = b.n_merge[clust] + 1;
            if b.is_basic[clust] {
                b.push_singleton(clust, gene, height[m - 1]);
            } else {
                on_branch[gene] = clust;
            }
            b.push_merging_height(clust, height[m - 1]);
            b.size[clust] += 1;
            b.n_merge[clust] = nm;
            b.n_singletons[clust] = ns;
            ind_merge_to_branch[m] = clust;
            root_branch = clust;
        } else {
            // two clusters
            let clust_a = ind_merge_to_branch[c1 as usize];
            let clust_b = ind_merge_to_branch[c2 as usize];
            let size_a = b.size[clust_a];
            let size_b = b.size[clust_b];
            // rank(sizes, ties.method="first") -> small = smaller (first on tie)
            let (mut small, mut large, sizes_small, sizes_large) = if size_a <= size_b {
                (clust_a, clust_b, size_a, size_b)
            } else {
                (clust_b, clust_a, size_b, size_a)
            };

            let small_ave_dist = if b.is_basic[small] {
                let cs = core_size(b.n_singletons[small], min_cluster_size);
                core_ave_dist(dist, &b.singletons[small], cs)
            } else {
                0.0
            };
            let large_ave_dist = if b.is_basic[large] {
                let cs = core_size(b.n_singletons[large], min_cluster_size);
                core_ave_dist(dist, &b.singletons[large], cs)
            } else {
                0.0
            };

            let smaller_scores = [
                b.is_basic[small],
                b.size[small] < min_cluster_size,
                small_ave_dist > max_abs_core_scatter,
                height[m - 1] - small_ave_dist < min_abs_gap,
                height[m - 1] < min_abs_split_height,
            ];
            let mut do_merge: bool;
            let mut smaller_fail_size: bool;
            if smaller_scores[0] && smaller_scores[1..].iter().any(|&x| x) {
                do_merge = true;
                smaller_fail_size = !(smaller_scores[2] || smaller_scores[3]);
            } else {
                let larger_scores = [
                    b.is_basic[large],
                    b.size[large] < min_cluster_size,
                    large_ave_dist > max_abs_core_scatter,
                    height[m - 1] - large_ave_dist < min_abs_gap,
                    height[m - 1] < min_abs_split_height,
                ];
                if larger_scores[0] && larger_scores[1..].iter().any(|&x| x) {
                    do_merge = true;
                    smaller_fail_size = !(larger_scores[2] || larger_scores[3]);
                    std::mem::swap(&mut small, &mut large);
                } else {
                    do_merge = false;
                    smaller_fail_size = false;
                }
            }

            if do_merge {
                b.fail_size[small] = smaller_fail_size;
                b.merged_into[small] = large;
                b.attach_height[small] = height[m - 1];
                b.is_top_basic[small] = false;
                let nss = b.n_singletons[small];
                let nsl = b.n_singletons[large];
                let ns = nss + nsl;
                if b.is_basic[large] {
                    let small_sings = b.singletons[small].clone();
                    for s in small_sings.iter().take(nss) {
                        b.push_singleton(large, *s, height[m - 1]);
                    }
                    b.n_singletons[large] = ns;
                } else {
                    if !b.is_basic[small] {
                        panic!("internal error: merging two composite clusters");
                    }
                    let small_sings = b.singletons[small].clone();
                    for s in small_sings.iter().take(nss) {
                        on_branch[*s] = large;
                    }
                }
                b.push_merging_height(large, height[m - 1]);
                b.size[large] = sizes_small + sizes_large;
                ind_merge_to_branch[m] = large;
                root_branch = large;
            } else {
                // not merging -> either split into a new composite branch or
                // fold into existing composite.
                if b.is_basic[large] && !b.is_basic[small] {
                    std::mem::swap(&mut small, &mut large);
                }
                if b.is_basic[large] || (pam_stage && pam_respects_dendro) {
                    n_branches += 1;
                    b.attach_height[large] = height[m - 1];
                    b.attach_height[small] = height[m - 1];
                    b.merged_into[large] = n_branches;
                    b.merged_into[small] = n_branches;
                    let mut add: Vec<usize> = Vec::new();
                    if b.is_basic[small] {
                        add.push(small);
                    } else {
                        add.extend_from_slice(&b.basic_clusters[small]);
                    }
                    if b.is_basic[large] {
                        add.push(large);
                    } else {
                        add.extend_from_slice(&b.basic_clusters[large]);
                    }
                    b.make_composite(n_branches, add, height[m - 1], sizes_small + sizes_large);
                    ind_merge_to_branch[m] = n_branches;
                    root_branch = n_branches;
                } else {
                    let mut add: Vec<usize> = if b.is_basic[small] {
                        vec![small]
                    } else {
                        b.basic_clusters[small].clone()
                    };
                    let nbl = b.n_basic_clusters[large];
                    let nb = nbl + add.len();
                    b.extend_basic_clusters(large, &mut add);
                    b.n_basic_clusters[large] = nb;
                    b.size[large] += b.size[small];
                    b.push_merging_height(large, height[m - 1]);
                    b.attach_height[small] = height[m - 1];
                    b.merged_into[small] = large;
                    ind_merge_to_branch[m] = large;
                    root_branch = large;
                }
            }
        }
    }

    // --- mark clusters ---
    let mut is_cluster = vec![false; n_branches + 1];
    let mut small_labels = vec![0usize; n_points + 1];
    for clust in 1..=n_branches {
        if b.attach_height[clust].is_nan() {
            b.attach_height[clust] = cut_height;
        }
        if b.is_top_basic[clust] {
            let coresize = core_size(b.n_singletons[clust], min_cluster_size);
            let core_scatter = core_ave_dist(dist, &b.singletons[clust], coresize);
            is_cluster[clust] = b.is_top_basic[clust]
                && b.size[clust] >= min_cluster_size
                && core_scatter < max_abs_core_scatter
                && b.attach_height[clust] - core_scatter > min_abs_gap;
        }
        if b.fail_size[clust] {
            for s in &b.singletons[clust] {
                small_labels[*s] = clust;
            }
        }
    }
    // respectSmallClusters = TRUE (default), so keep SmallLabels.

    // --- assign tree-cut labels ---
    let mut colors = vec![0i32; n_points + 1];
    let mut core_labels = vec![0i32; n_points + 1];
    let mut branch_labels = vec![0i32; n_branches + 1];
    let mut color = 0i32;
    for clust in 1..=n_branches {
        if is_cluster[clust] {
            color += 1;
            for s in &b.singletons[clust] {
                colors[*s] = color;
                small_labels[*s] = 0;
            }
            let coresize = core_size(b.n_singletons[clust], min_cluster_size);
            for s in b.singletons[clust].iter().take(coresize) {
                core_labels[*s] = color;
            }
            branch_labels[clust] = color;
        }
    }

    let labeled: Vec<usize> = (1..=n_points).filter(|&i| colors[i] != 0).collect();
    let mut unlabeled: Vec<usize> = (1..=n_points).filter(|&i| colors[i] == 0).collect();
    let unlabeled_exist = !unlabeled.is_empty();
    let n_proper_labels = if !labeled.is_empty() {
        unique_sorted(&colors.iter().map(|&i| i as usize).collect::<Vec<_>>(), |x| *x)
            .iter()
            .filter(|&&x| x != 0)
            .count() as i32
    } else {
        0
    };

    // --- PAM stage ---
    if pam_stage && unlabeled_exist && n_proper_labels > 0 {
        // useMedoids = FALSE branch.
        let mut cluster_diam = vec![0.0f64; (n_proper_labels + 1) as usize];
        for cluster in 1..=n_proper_labels {
            let in_cluster: Vec<usize> = (1..=n_points).filter(|&i| colors[i] == cluster).collect();
            let n_in_cluster = in_cluster.len();
            if n_in_cluster > 1 {
                let mut max_ave = 0.0f64;
                for &p in &in_cluster {
                    let mut s = 0.0;
                    for &q in &in_cluster {
                        s += dist[(p - 1, q - 1)];
                    }
                    let ave = s / (n_in_cluster as f64 - 1.0);
                    if ave > max_ave {
                        max_ave = ave;
                    }
                }
                cluster_diam[cluster as usize] = max_ave;
            } else {
                cluster_diam[cluster as usize] = 0.0;
            }
        }

        let colors_x = colors.clone();

        // Small clusters (respectSmallClusters = TRUE).
        if !small_labels.iter().any(|&s| s != 0) {
            // no small clusters; skip
        } else {
            let small_levs = unique_sorted(&small_labels.clone(), |x| *x);
            for &sclust in small_levs.iter().filter(|&&x| x != 0) {
                let in_cluster: Vec<usize> =
                    (1..=n_points).filter(|&i| small_labels[i] == sclust).collect();
                let on_br = unique_sorted(
                    &in_cluster.iter().map(|&i| on_branch[i]).collect::<Vec<_>>(),
                    |x| *x,
                );
                if on_br.len() > 1 {
                    panic!("internal error: small cluster on several branches");
                }
                let labels_on_branch: Vec<i32> = if on_br[0] > 0 {
                    let basic = &b.basic_clusters[on_br[0]];
                    basic.iter().map(|&bc| branch_labels[bc]).filter(|&l| l != 0).collect()
                } else {
                    Vec::new()
                };
                if labels_on_branch.is_empty() {
                    for &i in &in_cluster {
                        colors[i] = -1;
                    }
                    continue;
                }
                // mean distance from each small-cluster object to each label's objects
                let use_objects: Vec<usize> = (1..=n_points)
                    .filter(|&i| colors_x[i] != 0 && labels_on_branch.contains(&colors_x[i]))
                    .collect();
                if use_objects.is_empty() {
                    for &i in &in_cluster {
                        colors[i] = -1;
                    }
                    continue;
                }
                let mut mean_dist_by_label: Vec<(i32, f64)> = Vec::new();
                for &lab in labels_on_branch.iter() {
                    let objs: Vec<usize> = use_objects.iter().copied().filter(|&o| colors_x[o] == lab).collect();
                    if objs.is_empty() {
                        continue;
                    }
                    let mut s = 0.0;
                    for &o in &objs {
                        for &p in &in_cluster {
                            s += dist[(p - 1, o - 1)];
                        }
                    }
                    let mean = s / (objs.len() * in_cluster.len()) as f64;
                    mean_dist_by_label.push((lab, mean));
                }
                // nearest by mean
                let (nearest_label, nearest_dist) = mean_dist_by_label
                    .iter()
                    .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                    .map(|&(l, d)| (l, d))
                    .unwrap();
                if nearest_dist < cluster_diam[nearest_label as usize] || nearest_dist < max_pam_dist {
                    for &i in &in_cluster {
                        colors[i] = nearest_label;
                    }
                } else {
                    for &i in &in_cluster {
                        colors[i] = -1;
                    }
                }
            }
        }

        // Unlabeled objects (pamRespectsDendro = TRUE).
        unlabeled = (1..=n_points).filter(|&i| colors[i] == 0).collect();
        for &obj in unlabeled.iter().filter(|&&o| on_branch[o] > 0) {
            let on_br = on_branch[obj];
            let basic = b.basic_clusters[on_br].clone();
            let labels_on_branch: Vec<i32> =
                basic.iter().map(|&bc| branch_labels[bc]).filter(|&l| l != 0).collect();
            if labels_on_branch.is_empty() {
                continue;
            }
            let use_objects: Vec<usize> = (1..=n_points)
                .filter(|&i| colors_x[i] != 0 && labels_on_branch.contains(&colors_x[i]))
                .collect();
            if use_objects.is_empty() {
                continue;
            }
            let mut mean_dist_by_label: Vec<(i32, f64)> = Vec::new();
            for &lab in labels_on_branch.iter() {
                let objs: Vec<usize> = use_objects.iter().copied().filter(|&o| colors_x[o] == lab).collect();
                if objs.is_empty() {
                    continue;
                }
                let mut s = 0.0;
                for &o in &objs {
                    s += dist[(o - 1, obj - 1)];
                }
                let mean = s / objs.len() as f64;
                mean_dist_by_label.push((lab, mean));
            }
            let (nearest_label, nearest_dist) = mean_dist_by_label
                .iter()
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .map(|&(l, d)| (l, d))
                .unwrap();
            if nearest_dist < cluster_diam[nearest_label as usize] || nearest_dist < max_pam_dist {
                colors[obj] = nearest_label;
            }
        }
    }

    // --- renumber labels by size ---
    for c in colors.iter_mut() {
        if *c < 0 {
            *c = 0;
        }
    }
    let unlabeled_exist = colors[1..].iter().any(|&c| c == 0);

    // NumLabs = factor codes (1..k sorted), but we only need relative ordering.
    let uniq = unique_sorted(&colors[1..].to_vec(), |x| *x);
    // rank of cluster labels by size (descending), tie -> smaller label first
    let mut sizes_by_label: Vec<(i32, usize)> = Vec::new();
    for &u in &uniq {
        if u == 0 {
            continue;
        }
        let cnt = colors[1..].iter().filter(|&&c| c == u).count();
        sizes_by_label.push((u, cnt));
    }
    // sort by size desc, then label asc
    sizes_by_label.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let mut remap = std::collections::HashMap::new();
    for (new_label, &(old_label, _)) in sizes_by_label.iter().enumerate() {
        remap.insert(old_label, (new_label + 1) as i32);
    }

    let labels: Vec<i32> = (1..=n_points)
        .map(|i| {
            if colors[i] == 0 {
                0
            } else {
                *remap.get(&colors[i]).unwrap_or(&0)
            }
        })
        .collect();

    labels
}

fn core_ave_dist(dist: &Mat<f64>, core: &[usize], coresize: usize) -> f64 {
    if coresize <= 1 {
        return 0.0;
    }
    let c: Vec<usize> = core.iter().take(coresize).copied().collect();
    let mut s = 0.0;
    for &p in &c {
        let mut row = 0.0;
        for &q in &c {
            row += dist[(p - 1, q - 1)];
        }
        s += row / (coresize as f64 - 1.0);
    }
    s / coresize as f64
}

fn unique_sorted<T: Ord + Copy>(v: &[T], _f: fn(&T) -> T) -> Vec<T> {
    let mut u: Vec<T> = v.to_vec();
    u.sort();
    u.dedup();
    u
}

struct BranchState {
    singletons: Vec<Vec<usize>>,
    basic_clusters: Vec<Vec<usize>>,
    merging_heights: Vec<Vec<f64>>,
    size: Vec<usize>,
    n_merge: Vec<usize>,
    n_singletons: Vec<usize>,
    n_basic_clusters: Vec<usize>,
    is_basic: Vec<bool>,
    is_top_basic: Vec<bool>,
    fail_size: Vec<bool>,
    attach_height: Vec<f64>,
    merged_into: Vec<usize>,
}

impl BranchState {
    fn new(cap: usize) -> Self {
        let mut s = BranchState {
            singletons: vec![Vec::new()],
            basic_clusters: vec![Vec::new()],
            merging_heights: vec![Vec::new()],
            size: vec![0],
            n_merge: vec![0],
            n_singletons: vec![0],
            n_basic_clusters: vec![0],
            is_basic: vec![false],
            is_top_basic: vec![false],
            fail_size: vec![false],
            attach_height: vec![f64::NAN],
            merged_into: vec![0],
        };
        s.reserve(cap + 1);
        s
    }
    fn reserve(&mut self, cap: usize) {
        self.singletons.reserve(cap);
        self.basic_clusters.reserve(cap);
        self.merging_heights.reserve(cap);
        self.size.reserve(cap);
        self.n_merge.reserve(cap);
        self.n_singletons.reserve(cap);
        self.n_basic_clusters.reserve(cap);
        self.is_basic.reserve(cap);
        self.is_top_basic.reserve(cap);
        self.fail_size.reserve(cap);
        self.attach_height.reserve(cap);
        self.merged_into.reserve(cap);
    }
    fn set_basic(&mut self, idx: usize, g1: usize, g2: usize, h: f64) {
        if self.singletons.len() <= idx {
            self.push_empty();
        }
        self.singletons[idx] = vec![g1, g2];
        self.basic_clusters[idx] = Vec::new();
        self.merging_heights[idx] = vec![h, h];
        self.size[idx] = 2;
        self.n_merge[idx] = 1;
        self.n_singletons[idx] = 2;
        self.n_basic_clusters[idx] = 0;
        self.is_basic[idx] = true;
        self.is_top_basic[idx] = true;
        self.fail_size[idx] = false;
        self.attach_height[idx] = f64::NAN;
        self.merged_into[idx] = 0;
    }
    fn make_composite(&mut self, idx: usize, add: Vec<usize>, h: f64, size: usize) {
        if self.singletons.len() <= idx {
            self.push_empty();
        }
        self.singletons[idx] = Vec::new();
        self.basic_clusters[idx] = add.clone();
        self.merging_heights[idx] = vec![h, h];
        self.size[idx] = size;
        self.n_merge[idx] = 2;
        self.n_singletons[idx] = 0;
        self.n_basic_clusters[idx] = add.len();
        self.is_basic[idx] = false;
        self.is_top_basic[idx] = false;
        self.fail_size[idx] = false;
        self.attach_height[idx] = f64::NAN;
        self.merged_into[idx] = 0;
    }
    fn push_empty(&mut self) {
        self.singletons.push(Vec::new());
        self.basic_clusters.push(Vec::new());
        self.merging_heights.push(Vec::new());
        self.size.push(0);
        self.n_merge.push(0);
        self.n_singletons.push(0);
        self.n_basic_clusters.push(0);
        self.is_basic.push(false);
        self.is_top_basic.push(false);
        self.fail_size.push(false);
        self.attach_height.push(f64::NAN);
        self.merged_into.push(0);
    }
    fn push_singleton(&mut self, idx: usize, gene: usize, _h: f64) {
        self.singletons[idx].push(gene);
    }
    fn push_merging_height(&mut self, idx: usize, h: f64) {
        self.merging_heights[idx].push(h);
    }
    fn extend_basic_clusters(&mut self, idx: usize, add: &mut Vec<usize>) {
        self.basic_clusters[idx].append(add);
    }
}
