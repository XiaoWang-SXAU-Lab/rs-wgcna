use std::time::Instant;

use rs_wgcna::correlation::pearson;
use rs_wgcna::io::read_expression;
use rs_wgcna::network::{NetworkType, TomDenom, TomType};
use rs_wgcna::tom::tom_from_cor;

fn main() {
    let path = std::env::args().nth(1).expect("expr path");
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let exp = read_expression(&path).unwrap();
    eprintln!(
        "{} genes x {} samples, {} threads",
        exp.n_genes(),
        exp.n_samples(),
        threads
    );

    let t = Instant::now();
    let cor = pearson(&exp.data, threads).unwrap();
    let t_cor = t.elapsed();

    let t = Instant::now();
    let _tom = tom_from_cor(
        &cor,
        6.0,
        NetworkType::Unsigned,
        TomType::Signed,
        TomDenom::Min,
        threads,
    );
    let t_tom = t.elapsed();

    println!("cor_sec {:.4}", t_cor.as_secs_f64());
    println!("tom_sec {:.4}", t_tom.as_secs_f64());
}
