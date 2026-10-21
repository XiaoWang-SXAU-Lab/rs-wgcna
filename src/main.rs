//! Internal parity harness (not a user-facing CLI).
//!
//! Reads an expression matrix (plus optional colors / trait files) and writes
//! cor / adjacency / TOM / SFT / eigengenes / kME / module-trait as TSV so they
//! can be diffed against the R reference (scripts/r_reference.R).

use std::env;
use std::io::Write;

use faer::Mat;
use rs_wgcna::adjacency::adjacency;
use rs_wgcna::correlation::pearson;
use rs_wgcna::eigengene::module_eigengene;
use rs_wgcna::io::{read_expression, read_trait};
use rs_wgcna::kme::kme;
use rs_wgcna::module_trait::module_trait_cor;
use rs_wgcna::network::{NetworkType, TomDenom, TomType};
use rs_wgcna::sft::pick_soft_threshold;
use rs_wgcna::tom::tom_from_cor;

fn fmt(v: f64) -> String {
    if v.is_nan() {
        "NA".to_string()
    } else {
        format!("{v:.17}")
    }
}

fn write_mat(
    path: &str,
    m: &Mat<f64>,
    row_names: &[String],
    col_names: &[String],
) -> std::io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    write!(f, "ID")?;
    for c in col_names {
        write!(f, "\t{c}")?;
    }
    writeln!(f)?;
    for i in 0..m.nrows() {
        write!(f, "{}", row_names[i])?;
        for j in 0..m.ncols() {
            write!(f, "\t{}", fmt(m[(i, j)]))?;
        }
        writeln!(f)?;
    }
    Ok(())
}

fn read_lines(path: &str) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect()
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: rs-wgcna <expression.tsv> <power> <outdir> [colors.txt] [trait.tsv]");
        std::process::exit(2);
    }
    let expr_path = &args[1];
    let power: f64 = args[2].parse().expect("power");
    let outdir = &args[3];
    let colors_path = args.get(4);
    let trait_path = args.get(5);

    std::fs::create_dir_all(outdir).expect("create outdir");

    let exp = read_expression(expr_path).expect("read expression matrix");
    eprintln!(
        "loaded {} genes x {} samples",
        exp.n_genes(),
        exp.n_samples()
    );

    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let nt = NetworkType::Unsigned;

    // --- cor / adjacency / TOM ---
    let cor = pearson(&exp.data, threads).expect("correlation");
    write_mat(&format!("{outdir}/cor.tsv"), &cor, &exp.genes, &exp.genes).expect("write cor");

    let adj = adjacency(&cor, power, nt);
    write_mat(&format!("{outdir}/adjacency.tsv"), &adj, &exp.genes, &exp.genes)
        .expect("write adjacency");

    let tom = tom_from_cor(&cor, power, nt, TomType::Signed, TomDenom::Min, threads);
    write_mat(&format!("{outdir}/tom.tsv"), &tom, &exp.genes, &exp.genes).expect("write tom");

    // --- soft threshold (fit indices) ---
    let power_vector: Vec<f64> = (1..=10)
        .map(|p| p as f64)
        .chain((12..=20).step_by(2).map(|p| p as f64))
        .collect();
    let (_, sft_rows) = pick_soft_threshold(&cor, nt, &power_vector, 0.85, 10, false);
    {
        let mut f = std::fs::File::create(format!("{outdir}/sft.tsv")).expect("create sft");
        writeln!(f, "Power\tR.sq\tslope\ttruncated.R.sq\tmean.k\tmedian.k\tmax.k").unwrap();
        for r in &sft_rows {
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r[0], r[1], r[2], r[3], r[4], r[5], r[6]
            )
            .unwrap();
        }
    }

    // --- eigengenes + kME (require colors) ---
    if let Some(cp) = colors_path {
        let colors = read_lines(cp);
        assert_eq!(colors.len(), exp.n_genes(), "colors length != n_genes");
        let mut modules: Vec<&String> = colors.iter().collect();
        modules.sort();
        modules.dedup();

        let g = exp.n_genes();
        let n = exp.n_samples();
        let k = modules.len();
        let mut mes = Mat::<f64>::zeros(n, k);
        let mut mod_cols: Vec<String> = Vec::with_capacity(k);
        for (mi, m) in modules.iter().enumerate() {
            let idx: Vec<usize> = (0..g).filter(|&i| &colors[i] == *m).collect();
            let me = module_eigengene(&exp.data, &idx);
            for j in 0..n {
                mes[(j, mi)] = me[j];
            }
            mod_cols.push((*m).clone());
        }
        write_mat(
            &format!("{outdir}/eigengenes.tsv"),
            &mes,
            &exp.samples,
            &mod_cols,
        )
        .expect("write eigengenes");

        let km = kme(&exp.data, &mes);
        write_mat(&format!("{outdir}/kme.tsv"), &km, &exp.genes, &mod_cols).expect("write kme");

        // --- module-trait (require trait) ---
        if let Some(tp) = trait_path {
            let (_samples, traits, tdata) = read_trait(tp).expect("read trait");
            let (cor, p) = module_trait_cor(&mes, &tdata);
            write_mat(&format!("{outdir}/module_trait_cor.tsv"), &cor, &mod_cols, &traits)
                .expect("write mt cor");
            write_mat(&format!("{outdir}/module_trait_p.tsv"), &p, &mod_cols, &traits)
                .expect("write mt p");
        }
    }

    // --- full module detection (blockwiseModules single block) ---
    let min_module_size = 20usize.min(exp.n_genes() / 2);
    let module_colors = rs_wgcna::pipeline::blockwise_modules_single(
        &exp.data,
        power,
        nt,
        TomType::Signed,
        TomDenom::Min,
        min_module_size,
        0.15,
        2,
        0.995,
        threads,
    );
    {
        let mut f =
            std::fs::File::create(format!("{outdir}/module_colors.tsv")).expect("create colors");
        for c in &module_colors {
            writeln!(f, "{c}").unwrap();
        }
    }

    eprintln!("done.");
}
