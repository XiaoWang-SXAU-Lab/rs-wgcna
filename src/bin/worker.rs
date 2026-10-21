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
use rs_wgcna::filter::{
    low_expression_mask, normalize, top_genes_by, FilterMethod, NormalizeMethod,
};
use rs_wgcna::hclust::HClust;
use rs_wgcna::io::{read_expression, write_expression};
use rs_wgcna::network::{NetworkType, TomDenom, TomType};
use rs_wgcna::pipeline::blockwise_modules_with_tree;
use rs_wgcna::sft::pick_soft_threshold;

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
}

#[derive(Debug, Deserialize)]
#[serde(default)]
struct FilterParams {
    normalize_method: String, // raw | log10p1 | log
    n_keep: usize,
    filter_method: String, // mad | var
    sample_percentage: f64,
    cutoff: f64,
}

impl Default for FilterParams {
    fn default() -> Self {
        FilterParams {
            normalize_method: "log10p1".into(),
            n_keep: 8000,
            filter_method: "mad".into(),
            sample_percentage: 0.0,
            cutoff: 1.0,
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

// ---------------------------------------------------------------------------
// Output types
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct SftOutput {
    power_estimate: Option<f64>,
    network_type: String,
    fit_indices: Vec<FitIndex>,
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

fn read_power(job_dir: &Path, params: &NetworkParams) -> f64 {
    // Prefer the power recorded in step2_sft.json if present.
    if let Ok(s) = std::fs::read_to_string(job_dir.join("step2_sft.json")) {
        if let Ok(v) = serde_json::from_str::<serde_json::Value>(&s) {
            if let Some(p) = v.get("power_estimate").and_then(|x| x.as_f64()) {
                return p;
            }
        }
    }
    params.power
}

// ---------------------------------------------------------------------------
// Steps
// ---------------------------------------------------------------------------

fn run_filter(job_dir: &Path, params: &FilterParams) {
    write_progress(job_dir, "filter", "running", 0.0, "normalizing...");
    let exp = read_expression(job_dir.join("expr.tsv").to_str().unwrap()).expect("read expr.tsv");
    let norm = normalize(&exp.data, parse_norm(&params.normalize_method));

    write_progress(job_dir, "filter", "running", 0.4, "filtering low-expression genes...");
    let keep: Vec<bool> = if params.sample_percentage > 0.0 {
        low_expression_mask(&norm, params.cutoff, params.sample_percentage)
    } else {
        vec![true; exp.n_genes()]
    };

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

    write_progress(
        job_dir,
        "filter",
        "done",
        1.0,
        &format!("kept {} genes", filtered_genes.len()),
    );
}

fn run_soft_threshold(job_dir: &Path, params: &SoftThresholdParams, threads: usize) {
    write_progress(job_dir, "soft_threshold", "running", 0.2, "computing correlation...");
    let exp = read_filtered(job_dir);
    let nt = parse_network_type(&params.network_type);
    let cor = pearson(&exp.data, threads).expect("correlation");

    write_progress(job_dir, "soft_threshold", "running", 0.7, "scanning powers...");
    let power_vector: Vec<f64> = (1..=10)
        .map(|p| p as f64)
        .chain((12..=20).step_by(2).map(|p| p as f64))
        .collect();
    let (power_estimate, rows) =
        pick_soft_threshold(&cor, nt, &power_vector, params.rsquared_cut, 10, false);

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

    let out = SftOutput {
        power_estimate,
        network_type: params.network_type.clone(),
        fit_indices,
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
    let power = read_power(job_dir, params);

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
        eprintln!("usage: rs-wgcna-worker --job-dir <dir> --step <filter|soft_threshold|network> [--threads N]");
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
