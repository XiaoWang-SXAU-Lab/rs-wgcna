//! End-to-end single-block module detection, mirroring WGCNA `blockwiseModules`
//! (no block splitting, no missing-data filtering, no external branch splits).

use faer::Mat;

use crate::correlation::pearson;
use crate::cutree::cutree_hybrid;
use crate::eigengene::module_eigengene;
use crate::hclust::{hclust_average, HClust};
use crate::merge::{merge_close_modules, STANDARD_COLORS};
use crate::network::{NetworkType, TomDenom, TomType};
use crate::tom::tom_from_cor;

/// Convert integer labels to WGCNA colors (`labels2colors`, zeroIsGrey = TRUE).
pub fn labels_to_colors(labels: &[i32]) -> Vec<String> {
    labels
        .iter()
        .map(|&l| {
            if l == 0 {
                "grey".to_string()
            } else if (l as usize) <= STANDARD_COLORS.len() {
                STANDARD_COLORS[(l - 1) as usize].to_string()
            } else {
                format!("color{}", l)
            }
        })
        .collect()
}

/// Single-block module detection. Returns (final module colors, dendrogram).
#[allow(clippy::too_many_arguments)]
pub fn blockwise_modules_with_tree(
    expr: &Mat<f64>,
    power: f64,
    network_type: NetworkType,
    tom_type: TomType,
    tom_denom: TomDenom,
    min_module_size: usize,
    merge_cut_height: f64,
    deep_split: usize,
    detect_cut_height: f64,
    threads: usize,
) -> (Vec<String>, HClust) {
    let g = expr.nrows();

    // 1. correlation -> TOM -> dissimilarity.
    let cor = pearson(expr, threads).expect("correlation");
    let tom = tom_from_cor(&cor, power, network_type, tom_type, tom_denom, threads);
    let mut diss = tom.clone();
    for i in 0..g {
        for j in 0..g {
            diss[(i, j)] = 1.0 - tom[(i, j)];
        }
    }

    // 2. hierarchical clustering + dynamic tree cut.
    let tree = hclust_average(&diss);
    let labels = cutree_hybrid(
        &tree.merge,
        &tree.height,
        &diss,
        min_module_size,
        deep_split,
        detect_cut_height,
        true,
        true,
    );

    // 3. module eigengenes + kME-based cleanup (delete/trim modules).
    let mes = eigengenes_for_labels(expr, &labels);
    let mut labels = labels;
    cleanup_with_eigengenes(expr, &mut labels, min_module_size, network_type, &mes);

    // 3b. gene reassignment (Fisher p-value test, reassignThreshold = 1e-6).
    reassign_by_kme(expr, &mut labels, network_type, &mes, min_module_size);

    // 4. merge close modules + relabel.
    let colors = labels_to_colors(&labels);
    let merged = merge_close_modules(expr, &colors, merge_cut_height, threads);
    (merged, tree)
}

/// Single-block module detection. Returns final module colors (one per gene).
#[allow(clippy::too_many_arguments)]
pub fn blockwise_modules_single(
    expr: &Mat<f64>,
    power: f64,
    network_type: NetworkType,
    tom_type: TomType,
    tom_denom: TomDenom,
    min_module_size: usize,
    merge_cut_height: f64,
    deep_split: usize,
    detect_cut_height: f64,
    threads: usize,
) -> Vec<String> {
    blockwise_modules_with_tree(
        expr,
        power,
        network_type,
        tom_type,
        tom_denom,
        min_module_size,
        merge_cut_height,
        deep_split,
        detect_cut_height,
        threads,
    )
    .0
}

/// Remove modules with too few high-kME genes, and trim low-kME genes
/// (faithful to the `..checking kME in modules..` block in blockwiseModules).
pub fn cleanup_by_kme(
    expr: &Mat<f64>,
    labels: &mut [i32],
    min_module_size: usize,
    network_type: NetworkType,
) {
    let mes = eigengenes_for_labels(expr, labels);
    cleanup_with_eigengenes(expr, labels, min_module_size, network_type, &mes);
}

/// Columns retain the original numeric label mapping, including empty slots.
fn eigengenes_for_labels(expr: &Mat<f64>, labels: &[i32]) -> Mat<f64> {
    let n = expr.ncols();
    let max_label = *labels.iter().max().unwrap_or(&0) as usize;
    let mut mes = Mat::<f64>::zeros(n, max_label);
    for m in 1..=max_label {
        let idx: Vec<usize> = labels
            .iter()
            .enumerate()
            .filter(|(_, label)| **label == m as i32)
            .map(|(i, _)| i)
            .collect();
        if idx.is_empty() {
            continue;
        }
        let me = module_eigengene(expr, &idx);
        for j in 0..n {
            mes[(j, m - 1)] = me[j];
        }
    }
    mes
}

fn cleanup_with_eigengenes(
    expr: &Mat<f64>,
    labels: &mut [i32],
    min_module_size: usize,
    network_type: NetworkType,
    mes: &Mat<f64>,
) {
    let g = expr.nrows();
    let n = expr.ncols();
    let max_label = mes.ncols();
    let min_core_kme = 0.5;
    let min_core_kme_size = min_module_size as f64 / 3.0;
    let min_kme_to_stay = 0.3;

    let mut delete = vec![false; max_label + 1];
    for m in 1..=max_label {
        let mod_genes: Vec<usize> = (0..g).filter(|&i| labels[i] == m as i32).collect();
        if mod_genes.is_empty() {
            continue;
        }
        // kME = cor(gene, ME); unsigned network -> abs
        let kme: Vec<f64> = mod_genes
            .iter()
            .map(|&i| {
                let gene: Vec<f64> = (0..n).map(|j| expr[(i, j)]).collect();
                let me: Vec<f64> = (0..n).map(|j| mes[(j, m - 1)]).collect();
                let r = crate::stats::pearson_vec(&gene, &me);
                if network_type == NetworkType::Unsigned {
                    r.abs()
                } else {
                    r
                }
            })
            .collect();

        let n_high = kme.iter().filter(|k| **k > min_core_kme).count();
        if (n_high as f64) < min_core_kme_size {
            delete[m] = true;
        } else if kme.iter().any(|k| *k < min_kme_to_stay) {
            for (gi, k) in mod_genes.iter().zip(kme.iter()) {
                if *k < min_kme_to_stay {
                    labels[*gi] = 0;
                }
            }
            let remaining = mod_genes
                .iter()
                .zip(kme.iter())
                .filter(|(_, k)| **k >= min_kme_to_stay)
                .count();
            if remaining < min_module_size {
                delete[m] = true;
            }
        }
    }

    for i in 0..g {
        if labels[i] > 0 && delete[labels[i] as usize] {
            labels[i] = 0;
        }
    }
}

/// Reassign genes to a different module when their kME is significantly higher
/// there (Fisher p-value ratio < 1e-6), mirroring blockwiseModules'
/// "Check whether any of the already assigned genes should be re-assigned".
fn reassign_by_kme(
    expr: &Mat<f64>,
    labels: &mut [i32],
    network_type: NetworkType,
    mes: &Mat<f64>,
    min_module_size: usize,
) {
    use crate::stats::{cor_pvalue_fisher, pearson_vec};

    let g = expr.nrows();
    let n = expr.ncols();
    // WGCNA reuses the original dynamic-module eigengenes after cleanup;
    // deleted modules are excluded, but trimmed modules are not recomputed.
    let mut modules: Vec<i32> = labels.iter().copied().filter(|m| *m > 0).collect();
    modules.sort_unstable();
    modules.dedup();
    if modules.len() <= 1 {
        return;
    }
    let mut delete_modules = Vec::new();
    let reassign_threshold = 1e-6;
    let mut reassignments: Vec<(usize, i32)> = Vec::new();

    for &m in &modules {
        let before = reassignments.len();
        let mod_genes: Vec<usize> = (0..g).filter(|&i| labels[i] == m).collect();
        if mod_genes.is_empty() {
            continue;
        }
        for &gi in &mod_genes {
            let gene: Vec<f64> = (0..n).map(|j| expr[(gi, j)]).collect();
            let mut kme_module = f64::NAN;
            let mut best_col = 0usize;
            let mut kme_best = f64::NEG_INFINITY;
            for &module in &modules {
                let mm = (module - 1) as usize;
                let me: Vec<f64> = (0..n).map(|j| mes[(j, mm)]).collect();
                let r = pearson_vec(&gene, &me);
                let k = if network_type == NetworkType::Unsigned {
                    r.abs()
                } else {
                    r
                };
                if mm == (m - 1) as usize {
                    kme_module = k;
                }
                if k > kme_best {
                    kme_best = k;
                    best_col = module as usize;
                }
            }
            if kme_module.is_finite() && kme_module < kme_best {
                let p_module = cor_pvalue_fisher(kme_module, n);
                let p_best = cor_pvalue_fisher(kme_best, n);
                if p_best.is_finite() && p_module.is_finite() && p_best / p_module < reassign_threshold
                {
                    reassignments.push((gi, best_col as i32));
                }
            }
        }
        let moved = reassignments.len() - before;
        if moved > 0 && mod_genes.len() - moved < min_module_size {
            delete_modules.push(m);
        }
    }

    for (gi, new_label) in reassignments {
        labels[gi] = new_label;
    }
    for label in labels.iter_mut() {
        if delete_modules.contains(label) {
            *label = 0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reassignment_uses_supplied_dynamic_eigengenes() {
        // Centered orthogonal contrasts: R's supplied original eigengenes
        // give kME=(0.5,0.8), whose Fisher ratio at n=56 is below 1e-6.
        // Recomputing MEs from the singleton surviving modules would instead
        // give own kME=1 and incorrectly retain the first gene.
        let basis = |j: usize, bit: usize| if (j >> bit) & 1 == 0 { -1.0 } else { 1.0 };
        let mes = Mat::from_fn(56, 2, |j,m| basis(j,m));
        let expr = Mat::from_fn(2, 56, |i,j| if i == 0 {
            0.5*basis(j,0) + 0.8*basis(j,1) + 0.11_f64.sqrt()*basis(j,2)
        } else { basis(j,1) });
        let mut labels = vec![1,2];
        reassign_by_kme(&expr, &mut labels, NetworkType::Unsigned, &mes, 1);
        assert_eq!(labels, vec![2,2]);
    }

    #[test]
    fn reassignment_excludes_deleted_module_eigengenes() {
        let basis = |j: usize, bit: usize| if (j >> bit) & 1 == 0 { -1.0 } else { 1.0 };
        // Label 2 no longer exists after cleanup, despite its strong original ME.
        let mes = Mat::from_fn(56, 3, |j,m| basis(j,m));
        let expr = Mat::from_fn(2, 56, |i,j| if i == 0 {
            0.5*basis(j,0) + 0.8*basis(j,1) + 0.11_f64.sqrt()*basis(j,2)
        } else { basis(j,2) });
        let mut labels = vec![1,3];
        reassign_by_kme(&expr, &mut labels, NetworkType::Unsigned, &mes, 1);
        assert_eq!(labels, vec![1,3]);
    }
}
