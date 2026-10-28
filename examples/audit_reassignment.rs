//! TEST-ONLY kME threshold measurements on saved module-stage labels.
use faer::Mat;
use rs_wgcna::{eigengene::module_eigengene,filter::{normalize,NormalizeMethod},io::read_expression,stats::{pearson_vec,cor_pvalue_fisher}};
use std::collections::HashMap;
fn main(){
 let a:Vec<String>=std::env::args().collect();let exp=read_expression(&a[1]).unwrap();let norm=normalize(&exp.data,NormalizeMethod::Log10p1);
 let genes:Vec<String>=std::fs::read_to_string(format!("{}/filtered_genes.txt",a[2])).unwrap().lines().map(String::from).collect();
 let index:HashMap<&String,usize>=exp.genes.iter().enumerate().map(|(i,g)|(g,i)).collect();
 let expr=Mat::from_fn(genes.len(),norm.ncols(),|i,j|norm[(*index.get(&genes[i]).unwrap(),j)]);
 let gi=genes.iter().position(|g|g=="ENSG00000149541").unwrap();let gene:Vec<f64>=(0..expr.ncols()).map(|j|expr[(gi,j)]).collect();
 println!("stage\tgene\town_label\tbest_label\town_kME\tbest_kME\tp_own\tp_best\tratio\tthreshold\tdistance\treassign");
 for stage in ["dynamic","cleanup"] {
  let mut rdr=csv::ReaderBuilder::new().delimiter(b'\t').from_path(format!("{}/{stage}.tsv",a[2])).unwrap();
  let labels:Vec<i32>=rdr.records().enumerate().map(|(i,r)|{let r=r.unwrap();assert_eq!(&r[0],genes[i]);{let z:f64=r[1].parse().unwrap();assert_eq!(z.fract(),0.0);z as i32}}).collect();
  let mut mods:Vec<i32>=labels.iter().copied().filter(|&x|x>0).collect();mods.sort();mods.dedup();
  let k:Vec<f64>=mods.iter().map(|&m|{let idx:Vec<usize>=labels.iter().enumerate().filter(|(_,v)|**v==m).map(|(i,_)|i).collect();let me=module_eigengene(&expr,&idx);pearson_vec(&gene,&me).abs()}).collect();
  let own=mods.iter().position(|&m|m==labels[gi]).unwrap();let best=(0..k.len()).max_by(|&i,&j|k[i].partial_cmp(&k[j]).unwrap()).unwrap();let pown=cor_pvalue_fisher(k[own],expr.ncols());let pbest=cor_pvalue_fisher(k[best],expr.ncols());let ratio=pbest/pown;
  println!("{}\t{}\t{}\t{}\t{:.17}\t{:.17}\t{:.17e}\t{:.17e}\t{:.17e}\t1e-6\t{:.17e}\t{}",stage,genes[gi],labels[gi],mods[best],k[own],k[best],pown,pbest,ratio,ratio-1e-6,ratio<1e-6);
 }
}
