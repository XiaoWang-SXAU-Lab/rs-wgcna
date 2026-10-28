// TEST-ONLY diagnostic extension of examples/real_bench.rs.
//! Real-data benchmark: FPKM -> log10(x+1) -> MAD top-8000 -> WGCNA ->
//! module-trait. Writes results + per-step timings for the Quarto report.

use std::time::Instant;

use faer::Mat;
use rs_wgcna::correlation::pearson;
use rs_wgcna::eigengene::module_eigengene;
use rs_wgcna::filter::{normalize, top_genes_by, FilterMethod, NormalizeMethod};
use rs_wgcna::io::{read_expression, read_trait_aligned};
use rs_wgcna::module_trait::module_trait_cor;
use rs_wgcna::network::{NetworkType, TomDenom, TomType};
use rs_wgcna::pipeline::blockwise_modules_single;
use rs_wgcna::tom::tom_from_cor;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let fpkmpath = &args[1];
    let traitpath = &args[2];
    let outdir = &args[3];
    let n_keep: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(8000);
    std::fs::create_dir_all(outdir).unwrap();

    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let mut timings: Vec<(String, f64)> = Vec::new();

    // 1. read + normalize
    let t = Instant::now();
    let exp = read_expression(fpkmpath).unwrap();
    let norm = normalize(&exp.data, NormalizeMethod::Log10p1);
    timings.push(("log10_normalize".into(), t.elapsed().as_secs_f64()));

    // 2. MAD top-N
    let t = Instant::now();
    let idx = top_genes_by(&norm, n_keep, FilterMethod::Mad);
    let g = idx.len();
    let n = norm.ncols();
    let mut filtered = Mat::<f64>::zeros(g, n);
    for (r, &i) in idx.iter().enumerate() {
        for j in 0..n {
            filtered[(r, j)] = norm[(i, j)];
        }
    }
    let filtered_genes: Vec<String> = idx.iter().map(|&i| exp.genes[i].clone()).collect();
    timings.push(("mad_topN".into(), t.elapsed().as_secs_f64()));

    // 3. correlation
    let t = Instant::now();
    let cor = pearson(&filtered, threads).unwrap();
    timings.push(("cor".into(), t.elapsed().as_secs_f64()));

    // 4. TOM
    let t = Instant::now();
    let tom = tom_from_cor(
        &cor,
        6.0,
        NetworkType::Unsigned,
        TomType::Signed,
        TomDenom::Min,
        threads,
    );
    timings.push(("tom".into(), t.elapsed().as_secs_f64()));

    let probe: Vec<usize> = (0..g).step_by((g/64).max(1)).collect();
    let probe_genes: Vec<String> = probe.iter().map(|&i|filtered_genes[i].clone()).collect();
    for (name,mat) in [("cor_probe",&cor),("tom_probe",&tom)] {
        let sample = Mat::from_fn(probe.len(),g,|i,j|mat[(probe[i],j)]);
        write_mat(&format!("{outdir}/{name}.tsv"),&sample,&probe_genes,&filtered_genes);
    }
    let diss = Mat::from_fn(g,g,|i,j|1.0-tom[(i,j)]);
    let tree = rs_wgcna::hclust::hclust_average(&diss);
    let labels = rs_wgcna::cutree::cutree_hybrid(&tree.merge,&tree.height,&diss,30,2,0.995,true,true);
    let mut f = std::fs::File::create(format!("{outdir}/tree.tsv")).unwrap();
    writeln!(f,"ID\tleft\tright\theight").unwrap();
    for (i,((a,b),h)) in tree.merge.iter().zip(tree.height.iter()).enumerate() {writeln!(f,"{}\t{}\t{}\t{:.17}",i+1,a,b,h).unwrap();}
    let dynamic = Mat::from_fn(g,1,|i,_|labels[i] as f64);
    write_mat(&format!("{outdir}/dynamic.tsv"),&dynamic,&filtered_genes,&["label".into()]);
    let mut cleaned = labels.clone();
    rs_wgcna::pipeline::cleanup_by_kme(&filtered,&mut cleaned,30,NetworkType::Unsigned);
    let cleaned = Mat::from_fn(g,1,|i,_|cleaned[i] as f64);
    write_mat(&format!("{outdir}/cleanup.tsv"),&cleaned,&filtered_genes,&["label".into()]);
    let powers:Vec<f64>=(1..=10).map(|x|x as f64).chain((12..=20).step_by(2).map(|x|x as f64)).collect();
    let (_,sft)=rs_wgcna::sft::pick_soft_threshold(&cor,NetworkType::Unsigned,&powers,0.85,10,false);
    let sm=Mat::from_fn(sft.len(),6,|i,j|sft[i][j+1]);
    write_mat(&format!("{outdir}/sft.tsv"),&sm,&powers.iter().map(|p|p.to_string()).collect::<Vec<_>>(),&["R.sq","slope","truncated.R.sq","mean.k","median.k","max.k"].map(String::from));
    drop(diss);
    // 5. full module detection (cor+TOM+hclust+cutree+merge)
    let t = Instant::now();
    let colors = blockwise_modules_single(
        &filtered,
        6.0,
        NetworkType::Unsigned,
        TomType::Signed,
        TomDenom::Min,
        30,
        0.25,
        2,
        0.995,
        threads,
    );
    timings.push(("module_detection_total".into(), t.elapsed().as_secs_f64()));

    // 6. module eigengenes + module-trait
    let t = Instant::now();
    let (traits, tdata) = read_trait_aligned(traitpath, &exp.samples).expect("read aligned traits");
    // final module order (alphabetical non-grey)
    let mut mods: Vec<String> = colors.iter().cloned().collect();
    mods.sort();
    mods.dedup();
    let mods: Vec<String> = mods.into_iter().filter(|c| c != "grey").collect();
    let k = mods.len();
    let mut mes = Mat::<f64>::zeros(n, k);
    for (mi, m) in mods.iter().enumerate() {
        let midx: Vec<usize> = (0..g).filter(|&i| colors[i] == *m).collect();
        let me = module_eigengene(&filtered, &midx);
        for j in 0..n {
            mes[(j, mi)] = me[j];
        }
    }
    let (mt_cor, mt_p) = module_trait_cor(&mes, &tdata);
    timings.push(("module_trait".into(), t.elapsed().as_secs_f64()));

    write_mat(&format!("{outdir}/eigengenes.tsv"),&mes,&exp.samples,&mods);
    let km=rs_wgcna::kme::kme(&filtered,&mes);
    write_mat(&format!("{outdir}/kme.tsv"),&km,&filtered_genes,&mods);
    let intra=Mat::from_fn(g,1,|i,_|(0..g).filter(|&j|i!=j && colors[i]==colors[j]).map(|j|cor[(i,j)].abs().powi(6)).sum::<f64>());
    write_mat(&format!("{outdir}/intramodular.tsv"),&intra,&filtered_genes,&["within".into()]);
    // write outputs
    let mut f = std::fs::File::create(format!("{outdir}/filtered_genes.txt")).unwrap();
    for gi in &filtered_genes {
        use std::io::Write;
        writeln!(f, "{gi}").unwrap();
    }
    let mut f = std::fs::File::create(format!("{outdir}/module_colors.txt")).unwrap();
    for (gi, c) in filtered_genes.iter().zip(colors.iter()) {
        use std::io::Write;
        writeln!(f, "{gi}\t{c}").unwrap();
    }
    write_mat(
        &format!("{outdir}/module_trait_cor.tsv"),
        &mt_cor,
        &mods,
        &traits,
    );
    write_mat(
        &format!("{outdir}/module_trait_p.tsv"),
        &mt_p,
        &mods,
        &traits,
    );
    let mut f = std::fs::File::create(format!("{outdir}/timings.tsv")).unwrap();
    use std::io::Write;
    writeln!(f, "step\tseconds").unwrap();
    for (s, v) in &timings {
        writeln!(f, "{s}\t{v:.4}").unwrap();
    }
    let mut f = std::fs::File::create(format!("{outdir}/module_sizes.tsv")).unwrap();
    writeln!(f, "module\tsize").unwrap();
    for m in &mods {
        let sz = colors.iter().filter(|&c| c == m).count();
        writeln!(f, "{m}\t{sz}").unwrap();
    }

    eprintln!(
        "rust: {} genes -> {} MAD-top, {} modules, threads={}",
        exp.n_genes(),
        g,
        mods.len(),
        threads
    );
}

fn write_mat(path: &str, m: &Mat<f64>, rows: &[String], cols: &[String]) {
    use std::io::Write;
    let mut f = std::fs::File::create(path).unwrap();
    write!(f, "ID").unwrap();
    for c in cols {
        write!(f, "\t{c}").unwrap();
    }
    writeln!(f).unwrap();
    for i in 0..m.nrows() {
        write!(f, "{}", rows[i]).unwrap();
        for j in 0..m.ncols() {
            write!(f, "\t{:.17}", m[(i, j)]).unwrap();
        }
        writeln!(f).unwrap();
    }
}
