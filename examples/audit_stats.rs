//! TEST-ONLY independent statistical parameter grid.
fn main(){
 println!("ID\tr\tn\tstudent\tfisher");let mut i=0;
 for n in [3,4,5,10,20,40,56,100,1000] {
  for r in [0.0_f64,1e-8,1e-6,1e-4,0.001,0.005,0.01,0.02,0.05,0.1,0.25,0.5,0.9,0.99] {
   for sign in [-1.0,1.0] {i+=1;let x=r*sign;println!("{}\t{:.17}\t{}\t{:.17}\t{:.17}",i,x,n,rs_wgcna::stats::cor_pvalue(x,n),rs_wgcna::stats::cor_pvalue_fisher(x,n));}
  }
 }
}
