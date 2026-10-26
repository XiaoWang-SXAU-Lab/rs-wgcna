//! Internal worker: runs a single WGCNA pipeline step.
//!
//! Usage:
//!   rs-wgcna-worker --job-dir <dir> --step <filter|soft_threshold|network> [--threads N]
//!
//! The job directory holds the uploaded inputs and per-step outputs (which also
//! serve as the step-by-step cache). Progress is written to `<job-dir>/progress.json`.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use faer::Mat;
use serde::{Deserialize, Serialize};

use rs_wgcna::correlation::pearson;
use rs_wgcna::eigengene::module_eigengene;
use rs_wgcna::filter::{
    low_expression_mask, normalize, top_genes_by, FilterMethod, NormalizeMethod,
};
use rs_wgcna::hclust::{hclust_average, HClust};
use rs_wgcna::io::{read_expression, read_trait, write_expression};
use rs_wgcna::kme::kme;
use rs_wgcna::module_trait::module_trait_cor;
use rs_wgcna::network::{NetworkType, TomDenom, TomType};
use rs_wgcna::pipeline::blockwise_modules_with_tree;
use rs_wgcna::sft::{pick_soft_threshold, scale_free_check};
use rs_wgcna::stats::pearson_vec;
use rs_wgcna::tom::tom_from_cor;

// ---------------------------------------------------------------------------
// Parameters (read from <job-dir>/params.json)
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize, Default)]
struct Params {
    #[serde(default)]
    filter: FilterParams,
    #[serde(default)]
    soft_threshold: SoftThresholdParams,
    #[serde(default)]
    network: NetworkParams,
    #[serde(default)]
    module_trait: ModuleTraitParams,
    #[serde(default)]
    module_detail: ModuleDetailParams,
    #[serde(default)]
    hub: HubParams,
}

#[derive(Debug, Deserialize)]
#[serde(default)]
struct FilterParams {
    normalize_method: String, // raw | log10p1 | log
    n_keep: usize,
    filter_method: String, // mad | var
    sample_percentage: f64,
    cutoff: f64,
    exclude_samples: Vec<String>,
}

impl Default for FilterParams {
    fn default() -> Self {
        FilterParams {
            normalize_method: "log10p1".into(),
            n_keep: 8000,
            filter_method: "mad".into(),
            sample_percentage: 0.0,
            cutoff: 1.0,
            exclude_samples: Vec::new(),
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(default)]
struct SoftThresholdParams {
    network_type: String, // unsigned | signed | signed hybrid
    rsquared_cut: f64,
}

impl Default for SoftThresholdParams {
    fn default() -> Self {
        SoftThresholdParams {
            network_type: "unsigned".into(),
            rsquared_cut: 0.85,
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(default)]
struct NetworkParams {
    power: f64,
    network_type: String,
    tom_type: String,
    tom_denom: String,
    min_module_size: usize,
    merge_cut_height: f64,
    deep_split: usize,
    detect_cut_height: f64,
}

impl Default for NetworkParams {
    fn default() -> Self {
        NetworkParams {
            power: 6.0,
            network_type: "unsigned".into(),
            tom_type: "signed".into(),
            tom_denom: "min".into(),
            min_module_size: 30,
            merge_cut_height: 0.25,
            deep_split: 2,
            detect_cut_height: 0.995,
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(default)]
struct ModuleTraitParams {
    trait_source: String, // "uploaded"（上传的 trait 表）| "sample_info"（由 sample_info 生成）
    trait_column: String, // sample_info 中用于 one-hot 的列，如 "group" / "condition" / "time_h"
}

impl Default for ModuleTraitParams {
    fn default() -> Self {
        ModuleTraitParams {
            trait_source: "uploaded".into(),
            trait_column: String::new(),
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(default)]
struct ModuleDetailParams {
    select_module: String,
    select_trait: String,
}

impl Default for ModuleDetailParams {
    fn default() -> Self {
        ModuleDetailParams {
            select_module: String::new(),
            select_trait: String::new(),
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(default)]
struct HubParams {
    select_module: String,
    select_trait: String,
    kme_cutoff: f64,
    gs_cutoff: f64,
    weight_threshold: f64,
    power: f64,
    top_n_neighbors: usize,
    max_hubs: usize,
}

impl Default for HubParams {
    fn default() -> Self {
        HubParams {
            select_module: String::new(),
            select_trait: String::new(),
            kme_cutoff: 0.8,
            gs_cutoff: 0.0,
            weight_threshold: 0.02,
            power: 6.0,
            top_n_neighbors: 20,
            max_hubs: 30,
        }
    }
}

// ---------------------------------------------------------------------------
// Output types
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct SftOutput {
    power_estimate: Option<f64>,
    power_source: String, // "estimate"（R² 达标）| "experience"（经验兜底）
    network_type: String,
    fit_indices: Vec<FitIndex>,
    scale_free_check: ScaleFreeCheckOutput,
}

#[derive(Serialize)]
struct ScaleFreeCheckOutput {
    power: f64,
    hist_mids: Vec<f64>,
    hist_counts: Vec<f64>,
    log_k: Vec<f64>,
    log_p: Vec<f64>,
    intercept: f64,
    slope: f64,
    adj_r2: f64,
}

#[derive(Serialize)]
struct FitIndex {
    power: f64,
    rsq: f64,
    slope: f64,
    truncated_rsq: f64,
    mean_k: f64,
    median_k: f64,
    max_k: f64,
}

#[derive(Serialize)]
struct ModulesOutput {
    genes: Vec<String>,
    colors: Vec<String>,
    module_sizes: BTreeMap<String, usize>,
    dendrogram: DendrogramOutput,
}

#[derive(Serialize)]
struct DendrogramOutput {
    merge: Vec<(i32, i32)>,
    height: Vec<f64>,
    order: Vec<i32>,
}

impl From<&HClust> for DendrogramOutput {
    fn from(h: &HClust) -> Self {
        DendrogramOutput {
            merge: h.merge.clone(),
            height: h.height.clone(),
            order: h.order.clone(),
        }
    }
}

/// Minimal view of step3_modules.json needed by downstream steps.
#[derive(Deserialize, Default)]
struct ModulesJson {
    #[serde(default)]
    genes: Vec<String>,
    #[serde(default)]
    colors: Vec<String>,
}

#[derive(Serialize)]
struct MtaitOutput {
    modules: Vec<String>,
    traits: Vec<String>,
    cor: Vec<Vec<f64>>,
    p: Vec<Vec<f64>>,
    eigengenes: Vec<Vec<f64>>,
    genes: Vec<String>,
    kme: Vec<Vec<f64>>,
}

#[derive(Serialize)]
struct DetailOutput {
    module: String,
    trait_name: String,
    genes: Vec<String>,
    gs: Vec<f64>,
    mm: Vec<f64>,
    expression: Vec<Vec<f64>>,
    eigengene: Vec<f64>,
    samples: Vec<String>,
}

#[derive(Serialize)]
struct HubGene {
    gene: String,
    kme: f64,
    gs: f64,
}

#[derive(Serialize)]
struct NetNode {
    id: String,
    color: String,
    kme: f64,
    gs: f64,
    degree: usize,
    is_hub: bool,
}

#[derive(Serialize)]
struct NetEdge {
    source: String,
    target: String,
    weight: f64,
}

#[derive(Serialize)]
struct NetworkOutput {
    nodes: Vec<NetNode>,
    edges: Vec<NetEdge>,
}

#[derive(Serialize)]
struct HubOutput {
    module: String,
    trait_name: String,
    hub_genes: Vec<HubGene>,
    network: NetworkOutput,
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn write_progress(job_dir: &Path, step: &str, status: &str, progress: f64, message: &str) {
    let p = job_dir.join("progress.json");
    if let Ok(mut f) = std::fs::File::create(&p) {
        let json = serde_json::json!({
            "step": step,
            "status": status, // running | done | error
            "progress": progress,
            "message": message,
        });
        let _ = writeln!(f, "{json}");
    }
}

fn read_params(job_dir: &Path) -> Params {
    let path = job_dir.join("params.json");
    match std::fs::read_to_string(&path) {
        Ok(s) => serde_json::from_str(&s).unwrap_or_default(),
        Err(_) => Params::default(),
    }
}

fn parse_norm(s: &str) -> NormalizeMethod {
    match s {
        "raw" => NormalizeMethod::Raw,
        "log" => NormalizeMethod::Log,
        "log2p1" => NormalizeMethod::Log2p1,
        "log10" => NormalizeMethod::Log10,
        "log2" => NormalizeMethod::Log2,
        _ => NormalizeMethod::Log10p1,
    }
}

fn parse_filter_method(s: &str) -> FilterMethod {
    match s {
        "var" => FilterMethod::Var,
        _ => FilterMethod::Mad,
    }
}

fn parse_network_type(s: &str) -> NetworkType {
    NetworkType::from_str(s).unwrap_or(NetworkType::Unsigned)
}

fn parse_tom_type(s: &str) -> TomType {
    TomType::from_str(s).unwrap_or(TomType::Signed)
}

fn parse_tom_denom(s: &str) -> TomDenom {
    TomDenom::from_str(s).unwrap_or(TomDenom::Min)
}

fn read_filtered(job_dir: &Path) -> rs_wgcna::io::ExpData {
    let p = job_dir.join("step1_filtered.tsv");
    read_expression(p.to_str().unwrap()).expect("read step1_filtered.tsv")
}

// ---------------------------------------------------------------------------
// Steps
// ---------------------------------------------------------------------------

fn run_filter(job_dir: &Path, params: &FilterParams) {
    write_progress(job_dir, "filter", "running", 0.0, "normalizing...");
    let mut exp = read_expression(job_dir.join("expr.tsv").to_str().unwrap()).expect("read expr.tsv");
    let n_genes_before = exp.n_genes();
    let n_samples_before = exp.n_samples();

    // Remove excluded (outlier) samples before filtering.
    if !params.exclude_samples.is_empty() {
        let keep_cols: Vec<usize> = (0..exp.n_samples())
            .filter(|&j| !params.exclude_samples.iter().any(|s| *s == exp.samples[j]))
            .collect();
        let mut data = Mat::<f64>::zeros(exp.n_genes(), keep_cols.len());
        let mut samples = Vec::with_capacity(keep_cols.len());
        for (c, &j) in keep_cols.iter().enumerate() {
            for g in 0..exp.n_genes() {
                data[(g, c)] = exp.data[(g, j)];
            }
            samples.push(exp.samples[j].clone());
        }
        exp.data = data;
        exp.samples = samples;
    }

    // vst（DESeq2 方差稳定变换，纯 Rust）用于 raw count；其余用元素级归一化。
    let norm = if params.normalize_method == "vst" {
        rs_wgcna::vst::vst_transform(&exp.data)
    } else {
        normalize(&exp.data, parse_norm(&params.normalize_method))
    };

    write_progress(job_dir, "filter", "running", 0.4, "filtering low-expression genes...");
    let keep: Vec<bool> = if params.sample_percentage > 0.0 {
        low_expression_mask(&norm, params.cutoff, params.sample_percentage)
    } else {
        vec![true; exp.n_genes()]
    };
    let n_genes_after_low = keep.iter().filter(|&&k| k).count();

    write_progress(job_dir, "filter", "running", 0.7, "selecting top genes...");
    // Select top-N by MAD/Var among the kept genes.
    let (filtered, filtered_genes) = {
        let kept_idx: Vec<usize> = (0..exp.n_genes()).filter(|&i| keep[i]).collect();
        let g = kept_idx.len();
        let n = norm.ncols();
        let mut kept_mat = Mat::<f64>::zeros(g, n);
        for (r, &i) in kept_idx.iter().enumerate() {
            for j in 0..n {
                kept_mat[(r, j)] = norm[(i, j)];
            }
        }
        let kept_names: Vec<String> = kept_idx.iter().map(|&i| exp.genes[i].clone()).collect();
        let top = top_genes_by(&kept_mat, params.n_keep.min(g), parse_filter_method(&params.filter_method));
        let mut out = Mat::<f64>::zeros(top.len(), n);
        let mut out_names = Vec::with_capacity(top.len());
        for (r, &ti) in top.iter().enumerate() {
            for j in 0..n {
                out[(r, j)] = kept_mat[(ti, j)];
            }
            out_names.push(kept_names[ti].clone());
        }
        (out, out_names)
    };

    write_expression(
        job_dir.join("step1_filtered.tsv").to_str().unwrap(),
        &filtered,
        &filtered_genes,
        &exp.samples,
    )
    .expect("write step1_filtered.tsv");

    let mut f = std::fs::File::create(job_dir.join("step1_genes.txt")).expect("create genes.txt");
    for gi in &filtered_genes {
        writeln!(f, "{gi}").unwrap();
    }

    // Sample clustering tree (average-linkage on Euclidean sample distance).
    write_progress(job_dir, "filter", "running", 0.9, "clustering samples...");
    let ns = filtered.ncols();
    let mut diss = Mat::<f64>::zeros(ns, ns);
    for i in 0..ns {
        for j in 0..ns {
            let mut s = 0.0;
            for g in 0..filtered.nrows() {
                let d = filtered[(g, i)] - filtered[(g, j)];
                s += d * d;
            }
            diss[(i, j)] = s.sqrt();
        }
    }
    let tree = hclust_average(&diss);

    let out = serde_json::json!({
        "n_genes_before": n_genes_before,
        "n_genes_after_low": n_genes_after_low,
        "n_genes_final": filtered_genes.len(),
        "n_samples_before": n_samples_before,
        "n_samples": exp.samples.len(),
        "samples": exp.samples,
        "sample_tree": { "merge": tree.merge, "height": tree.height, "order": tree.order },
    });
    write_json(job_dir, "step1_filter.json", &out);

    write_progress(
        job_dir,
        "filter",
        "done",
        1.0,
        &format!("kept {} genes", filtered_genes.len()),
    );
}

/// WGCNA FAQ 经验 power 值：R² 未达标（powerEstimate 为 NA）时兜底。
/// 样本数越少，经验 power 越大（相关性估计噪声更大，需更强幂次）。
fn experience_power(n_samples: usize, nt: NetworkType) -> f64 {
    let signed = matches!(nt, NetworkType::Signed | NetworkType::SignedHybrid);
    if n_samples < 20 {
        if signed { 18.0 } else { 9.0 }
    } else if n_samples < 30 {
        if signed { 16.0 } else { 8.0 }
    } else if n_samples < 40 {
        if signed { 14.0 } else { 7.0 }
    } else {
        if signed { 12.0 } else { 6.0 }
    }
}

fn run_soft_threshold(job_dir: &Path, params: &SoftThresholdParams, threads: usize) {
    write_progress(job_dir, "soft_threshold", "running", 0.2, "computing correlation...");
    let exp = read_filtered(job_dir);
    let nt = parse_network_type(&params.network_type);
    let cor = pearson(&exp.data, threads).expect("correlation");

    write_progress(job_dir, "soft_threshold", "running", 0.7, "scanning powers...");
    // 扫满 1..=20 全部 power，不提前停止。
    let power_vector: Vec<f64> = (1..=20).map(|p| p as f64).collect();
    let (power_estimate, rows) =
        pick_soft_threshold(&cor, nt, &power_vector, params.rsquared_cut, 10, false);

    // R² 未达标时，用 WGCNA FAQ 经验 power 兜底（同步自 Shiny app 的 getpower 逻辑）。
    let n_samples = exp.data.ncols();
    let (power_estimate, power_source) = match power_estimate {
        Some(p) => (Some(p), "estimate".to_string()),
        None => {
            let ep = experience_power(n_samples, nt);
            (Some(ep), "experience".to_string())
        }
    };

    let fit_indices = rows
        .iter()
        .map(|r| FitIndex {
            power: r[0],
            rsq: r[1],
            slope: r[2],
            truncated_rsq: r[3],
            mean_k: r[4],
            median_k: r[5],
            max_k: r[6],
        })
        .collect();

    // Scale-free topology check for the selected power (connectivity histogram
    // + log-log regression), mirroring the old Shiny powertest.
    let selected_power = power_estimate.expect("power estimate");
    let sfc = scale_free_check(&cor, selected_power, nt, 10, false);

    let out = SftOutput {
        power_estimate,
        power_source,
        network_type: params.network_type.clone(),
        fit_indices,
        scale_free_check: ScaleFreeCheckOutput {
            power: selected_power,
            hist_mids: sfc.hist_mids,
            hist_counts: sfc.hist_counts,
            log_k: sfc.log_k,
            log_p: sfc.log_p,
            intercept: sfc.intercept,
            slope: sfc.slope,
            adj_r2: sfc.adj_r2,
        },
    };
    let s = serde_json::to_string_pretty(&out).unwrap();
    std::fs::write(job_dir.join("step2_sft.json"), s).expect("write step2_sft.json");

    write_progress(job_dir, "soft_threshold", "done", 1.0, "soft threshold selected");
}

fn run_network(job_dir: &Path, params: &NetworkParams, threads: usize) {
    write_progress(job_dir, "network", "running", 0.1, "computing TOM...");
    let exp = read_filtered(job_dir);
    let nt = parse_network_type(&params.network_type);
    let tom_type = parse_tom_type(&params.tom_type);
    let tom_denom = parse_tom_denom(&params.tom_denom);
    let power = params.power;

    write_progress(job_dir, "network", "running", 0.6, "clustering & cutting modules...");
    let (colors, tree) = blockwise_modules_with_tree(
        &exp.data,
        power,
        nt,
        tom_type,
        tom_denom,
        params.min_module_size,
        params.merge_cut_height,
        params.deep_split,
        params.detect_cut_height,
        threads,
    );

    let mut module_sizes = BTreeMap::new();
    for c in &colors {
        *module_sizes.entry(c.clone()).or_insert(0) += 1;
    }

    let out = ModulesOutput {
        genes: exp.genes.clone(),
        colors,
        module_sizes,
        dendrogram: DendrogramOutput::from(&tree),
    };
    let s = serde_json::to_string(&out).unwrap();
    std::fs::write(job_dir.join("step3_modules.json"), s).expect("write step3_modules.json");

    write_progress(job_dir, "network", "done", 1.0, "module detection done");
}

fn read_json<T: for<'de> Deserialize<'de>>(job_dir: &Path, name: &str) -> T {
    let s = std::fs::read_to_string(job_dir.join(name)).unwrap_or_else(|_| "{}".into());
    serde_json::from_str(&s).unwrap_or_else(|e| panic!("invalid {name}: {e}"))
}

fn write_json<T: Serialize>(job_dir: &Path, name: &str, v: &T) {
    let s = serde_json::to_string(v).expect("serialize");
    std::fs::write(job_dir.join(name), s).expect("write output");
}

fn read_modules(job_dir: &Path) -> ModulesJson {
    read_json(job_dir, "step3_modules.json")
}

fn sorted_modules(colors: &[String]) -> Vec<String> {
    let mut mods: Vec<String> = colors.to_vec();
    mods.sort();
    mods.dedup();
    mods.into_iter().filter(|c| c != "grey").collect()
}

/// Fallback module selection for one-step mode: the largest non-grey module.
fn largest_module(colors: &[String]) -> String {
    let mut sizes: BTreeMap<&str, usize> = BTreeMap::new();
    for c in colors {
        *sizes.entry(c.as_str()).or_insert(0) += 1;
    }
    sizes
        .into_iter()
        .filter(|(c, _)| *c != "grey")
        .max_by_key(|(_, n)| *n)
        .map(|(c, _)| c.to_string())
        .unwrap_or_else(|| "grey".to_string())
}

fn trait_column(traits: &[String], name: &str) -> usize {
    traits.iter().position(|t| t == name).unwrap_or(0)
}

/// Build a one-hot trait matrix from sample_info.tsv's `column`, aligned to
/// `exp_samples`. Returns (category_names, samples × categories one-hot matrix).
fn sample_info_one_hot(
    path: &str,
    column: &str,
    exp_samples: &[String],
) -> rs_wgcna::error::Result<(Vec<String>, Mat<f64>)> {
    use std::collections::HashMap;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_path(path)?;
    let headers = rdr.headers()?.clone();
    let col_idx = headers.iter().position(|h| h == column).ok_or_else(|| {
        rs_wgcna::error::WgcnaError::InvalidInput(format!("sample_info column '{column}' not found"))
    })?;

    let mut map: HashMap<String, String> = HashMap::new();
    for rec in rdr.records() {
        let rec = rec?;
        let id = rec.get(0).unwrap_or("").to_string();
        let val = rec.get(col_idx).unwrap_or("").to_string();
        map.insert(id, val);
    }

    let mut categories: Vec<String> = Vec::new();
    for s in exp_samples {
        let v = map.get(s).cloned().unwrap_or_default();
        if !categories.contains(&v) {
            categories.push(v);
        }
    }

    let n = exp_samples.len();
    let c = categories.len();
    let mut data = Mat::<f64>::zeros(n, c);
    for i in 0..n {
        let v = map.get(&exp_samples[i]).cloned().unwrap_or_default();
        let ci = categories.iter().position(|x| *x == v).unwrap_or(0);
        data[(i, ci)] = 1.0;
    }
    Ok((categories, data))
}

fn run_module_trait(job_dir: &Path, params: &ModuleTraitParams, _threads: usize) {
    write_progress(job_dir, "module_trait", "running", 0.2, "reading modules...");
    let exp = read_filtered(job_dir);
    let mj = read_modules(job_dir);
    // trait source: uploaded table OR generated from sample_info
    let (traits, tdata) = if params.trait_source == "sample_info" {
        sample_info_one_hot(
            job_dir.join("sample_info.tsv").to_str().unwrap(),
            &params.trait_column,
            &exp.samples,
        )
        .expect("read sample_info.tsv")
    } else {
        let (_samples, traits, tdata) =
            read_trait(job_dir.join("trait.tsv").to_str().unwrap()).expect("read trait.tsv");
        (traits, tdata)
    };

    let mods = sorted_modules(&mj.colors);
    let g = mj.genes.len();
    let n = exp.n_samples();
    let mut mes = Mat::<f64>::zeros(n, mods.len());
    for (mi, m) in mods.iter().enumerate() {
        let idx: Vec<usize> = (0..g).filter(|&i| mj.colors[i] == *m).collect();
        let me = module_eigengene(&exp.data, &idx);
        for j in 0..n {
            mes[(j, mi)] = me[j];
        }
    }

    write_progress(job_dir, "module_trait", "running", 0.7, "correlating with traits...");
    let (cor, p) = module_trait_cor(&mes, &tdata);
    let cor_v: Vec<Vec<f64>> = (0..cor.nrows())
        .map(|i| (0..cor.ncols()).map(|j| cor[(i, j)]).collect())
        .collect();
    let p_v: Vec<Vec<f64>> = (0..p.nrows())
        .map(|i| (0..p.ncols()).map(|j| p[(i, j)]).collect())
        .collect();
    let me_v: Vec<Vec<f64>> = (0..mes.nrows())
        .map(|i| (0..mes.ncols()).map(|j| mes[(i, j)]).collect())
        .collect();
    let kme_mat = kme(&exp.data, &mes);
    let kme_v: Vec<Vec<f64>> = (0..kme_mat.nrows())
        .map(|i| (0..kme_mat.ncols()).map(|j| kme_mat[(i, j)]).collect())
        .collect();

    write_json(
        job_dir,
        "step4_mtrait.json",
        &MtaitOutput {
            modules: mods,
            traits,
            cor: cor_v,
            p: p_v,
            eigengenes: me_v,
            genes: mj.genes.clone(),
            kme: kme_v,
        },
    );
    write_progress(job_dir, "module_trait", "done", 1.0, "module-trait done");
}

fn run_module_detail(job_dir: &Path, params: &ModuleDetailParams, _threads: usize) {
    let exp = read_filtered(job_dir);
    let mj = read_modules(job_dir);
    let (_s, traits, tdata) =
        read_trait(job_dir.join("trait.tsv").to_str().unwrap()).expect("read trait.tsv");
    let t_idx = trait_column(&traits, &params.select_trait);
    let trait_v: Vec<f64> = (0..tdata.nrows()).map(|j| tdata[(j, t_idx)]).collect();

    let sel_module = if params.select_module.is_empty() {
        largest_module(&mj.colors)
    } else {
        params.select_module.clone()
    };
    let g = mj.genes.len();
    let mod_idx: Vec<usize> = (0..g).filter(|&i| mj.colors[i] == sel_module).collect();
    let me = module_eigengene(&exp.data, &mod_idx);

    let mut genes = Vec::new();
    let mut gs = Vec::new();
    let mut mm = Vec::new();
    let mut expression = Vec::new();
    for &gi in &mod_idx {
        let gene: Vec<f64> = (0..exp.n_samples()).map(|j| exp.data[(gi, j)]).collect();
        genes.push(mj.genes[gi].clone());
        gs.push(pearson_vec(&gene, &trait_v));
        mm.push(pearson_vec(&gene, &me));
        expression.push(gene);
    }

    write_json(
        job_dir,
        "step5_detail.json",
        &DetailOutput {
            module: sel_module.clone(),
            trait_name: traits.get(t_idx).cloned().unwrap_or_default(),
            genes,
            gs,
            mm,
            expression,
            eigengene: me,
            samples: exp.samples.clone(),
        },
    );
    write_progress(job_dir, "module_detail", "done", 1.0, "module detail done");
}

fn run_hub(job_dir: &Path, params: &HubParams, net_params: &NetworkParams, threads: usize) {
    let exp = read_filtered(job_dir);
    let mj = read_modules(job_dir);
    let (_s, traits, tdata) =
        read_trait(job_dir.join("trait.tsv").to_str().unwrap()).expect("read trait.tsv");
    let t_idx = trait_column(&traits, &params.select_trait);
    let trait_v: Vec<f64> = (0..tdata.nrows()).map(|j| tdata[(j, t_idx)]).collect();

    let sel_module = if params.select_module.is_empty() {
        largest_module(&mj.colors)
    } else {
        params.select_module.clone()
    };
    let g = mj.genes.len();
    let mod_idx: Vec<usize> = (0..g).filter(|&i| mj.colors[i] == sel_module).collect();
    let me = module_eigengene(&exp.data, &mod_idx);

    use std::collections::{HashMap, HashSet};

    // Score every module gene by kME + GS once, then reuse for the table, the
    // hub seeds, and the candidate pool.
    let mut scored: Vec<(usize, f64, f64)> = mod_idx
        .iter()
        .map(|&gi| {
            let gene: Vec<f64> = (0..exp.n_samples()).map(|j| exp.data[(gi, j)]).collect();
            let kme = pearson_vec(&gene, &me).abs();
            let gs = pearson_vec(&gene, &trait_v).abs();
            (gi, kme, gs)
        })
        .collect();
    scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    let metric: HashMap<usize, (f64, f64)> = scored.iter().map(|&(gi, k, g)| (gi, (k, g))).collect();

    // Hub gene table: |kME| > cutoff AND |GS| > cutoff, sorted by kME (desc).
    let mut hub = Vec::new();
    for &(gi, k, g) in &scored {
        if k > params.kme_cutoff && g > params.gs_cutoff {
            hub.push(HubGene {
                gene: mj.genes[gi].clone(),
                kme: k,
                gs: g,
            });
        }
    }

    // --- network visualization: hub genes (kME > cutoff) + their neighbors ---

    // 2. hub genes = kME > cutoff, limited to max_hubs
    let hub_genes: Vec<usize> = scored
        .iter()
        .filter(|&&(_, k, _)| k > params.kme_cutoff)
        .take(params.max_hubs)
        .map(|&(gi, _, _)| gi)
        .collect();

    // 3. candidates = hubs + top pool by kME (for TOM computation)
    const CANDIDATE_POOL: usize = 400;
    let mut candidate: Vec<usize> = hub_genes.clone();
    let mut seen: HashSet<usize> = hub_genes.iter().copied().collect();
    for &(gi, _, _) in &scored {
        if candidate.len() >= CANDIDATE_POOL {
            break;
        }
        if seen.insert(gi) {
            candidate.push(gi);
        }
    }
    let ng = candidate.len();
    let gene_to_idx: HashMap<usize, usize> = candidate.iter().enumerate().map(|(i, &gi)| (gi, i)).collect();
    let hub_pos: HashSet<usize> = hub_genes.iter().map(|&gi| gene_to_idx[&gi]).collect();

    // 4. TOM over candidates
    let mut cor = Mat::<f64>::zeros(ng, ng);
    for a in 0..ng {
        let ga: Vec<f64> = (0..exp.n_samples()).map(|j| exp.data[(candidate[a], j)]).collect();
        for b in (a + 1)..ng {
            let gb: Vec<f64> = (0..exp.n_samples()).map(|j| exp.data[(candidate[b], j)]).collect();
            let c = pearson_vec(&ga, &gb);
            cor[(a, b)] = c;
            cor[(b, a)] = c;
        }
        cor[(a, a)] = 1.0;
    }
    let tom = tom_from_cor(
        &cor,
        params.power,
        parse_network_type(&net_params.network_type),
        parse_tom_type(&net_params.tom_type),
        parse_tom_denom(&net_params.tom_denom),
        threads,
    );

    // 5. per-hub top-N neighbors (by weight, above threshold)
    let mut kept: HashSet<(usize, usize)> = HashSet::new();
    for &h in &hub_pos {
        let mut nb: Vec<(usize, f64)> = Vec::new();
        for b in 0..ng {
            if b == h {
                continue;
            }
            let w = tom[(h, b)];
            if w > params.weight_threshold {
                nb.push((b, w));
            }
        }
        nb.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        for &(b, _) in nb.iter().take(params.top_n_neighbors) {
            kept.insert(if h < b { (h, b) } else { (b, h) });
        }
    }

    // 6. nodes (hubs + neighbors) + edges + degree
    let mut node_set: HashSet<usize> = hub_pos.clone();
    let mut edges: Vec<NetEdge> = Vec::new();
    let mut degree: HashMap<usize, usize> = HashMap::new();
    for &(a, b) in &kept {
        node_set.insert(a);
        node_set.insert(b);
        *degree.entry(a).or_insert(0) += 1;
        *degree.entry(b).or_insert(0) += 1;
        edges.push(NetEdge {
            source: mj.genes[candidate[a]].clone(),
            target: mj.genes[candidate[b]].clone(),
            weight: tom[(a, b)],
        });
    }

    let mut nodes: Vec<NetNode> = Vec::new();
    for &idx in &node_set {
        let gi = candidate[idx];
        let (k, g) = metric.get(&gi).copied().unwrap_or((0.0, 0.0));
        nodes.push(NetNode {
            id: mj.genes[gi].clone(),
            color: mj.colors[gi].clone(),
            kme: k,
            gs: g,
            degree: degree.get(&idx).copied().unwrap_or(0),
            is_hub: hub_pos.contains(&idx),
        });
    }

    write_json(
        job_dir,
        "step6_hub.json",
        &HubOutput {
            module: sel_module.clone(),
            trait_name: traits.get(t_idx).cloned().unwrap_or_default(),
            hub_genes: hub,
            network: NetworkOutput { nodes, edges },
        },
    );
    write_progress(job_dir, "hub", "done", 1.0, "hub genes done");
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut job_dir: Option<PathBuf> = None;
    let mut step: Option<String> = None;
    let mut threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--job-dir" => {
                job_dir = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--step" => {
                step = Some(args[i + 1].clone());
                i += 2;
            }
            "--threads" => {
                threads = args[i + 1].parse().unwrap_or(threads);
                i += 2;
            }
            _ => i += 1,
        }
    }

    let job_dir = job_dir.unwrap_or_else(|| {
        eprintln!("usage: rs-wgcna-worker --job-dir <dir> --step <filter|soft_threshold|network|module_trait|module_detail|hub> [--threads N]");
        std::process::exit(2);
    });
    let step = step.unwrap_or_else(|| {
        eprintln!("missing --step");
        std::process::exit(2);
    });

    let params = read_params(&job_dir);

    let result = std::panic::catch_unwind(|| match step.as_str() {
        "filter" => run_filter(&job_dir, &params.filter),
        "soft_threshold" => run_soft_threshold(&job_dir, &params.soft_threshold, threads),
        "network" => run_network(&job_dir, &params.network, threads),
        "module_trait" => run_module_trait(&job_dir, &params.module_trait, threads),
        "module_detail" => run_module_detail(&job_dir, &params.module_detail, threads),
        "hub" => run_hub(&job_dir, &params.hub, &params.network, threads),
        other => {
            eprintln!("unknown step: {other}");
            std::process::exit(2);
        }
    });

    if let Err(e) = result {
        let msg = if let Some(s) = e.downcast_ref::<&str>() {
            s.to_string()
        } else if let Some(s) = e.downcast_ref::<String>() {
            s.clone()
        } else {
            "worker panic".to_string()
        };
        write_progress(&job_dir, &step, "error", 1.0, &msg);
        eprintln!("worker error: {msg}");
        std::process::exit(1);
    }
}
