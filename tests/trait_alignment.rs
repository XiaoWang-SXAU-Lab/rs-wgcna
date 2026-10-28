use faer::Mat;
use rs_wgcna::io::{align_trait_samples, read_trait_aligned};
fn names(ids: &[&str])->Vec<String>{ids.iter().map(|s|s.to_string()).collect()}
#[test]
fn align_trait_rows_by_sample_id() {
    let data=Mat::from_fn(3,2,|i,j|(10*i+j) as f64);
    let aligned=align_trait_samples(&data,&names(&["c","a","b"]),&names(&["a","b","c"])).unwrap();
    assert_eq!((aligned[(0,0)],aligned[(1,1)],aligned[(2,0)]),(10.0,21.0,0.0));
}
#[test]
fn reject_ambiguous_or_different_sample_sets() {
    let data=Mat::zeros(2,1);
    for (source,target) in [(["a","a"],["a","b"]),(["a","b"],["a","a"]),(["a","x"],["a","b"]),(["a",""],["a","b"])] {
        assert!(align_trait_samples(&data,&names(&source),&names(&target)).is_err());
    }
    assert!(align_trait_samples(&data,&names(&["a","b"]),&names(&["a"])).is_err());
    assert!(align_trait_samples(&Mat::zeros(1,1),&names(&["a","b"]),&names(&["a","b"])).is_err());
}
#[test]
fn reject_missing_and_extra_trait_cells() {
    let base=std::env::temp_dir().join(format!("rswgcna-trait-regression-{}",std::process::id()));
    std::fs::create_dir_all(&base).unwrap();
    for (i,text) in ["ID\tt1\tt2\na\t1\n", "ID\tt1\na\t1\t2\n", "ID\tt1\na\t1\na\t2\n"].iter().enumerate() {
        let p=base.join(format!("{i}.tsv"));std::fs::write(&p,text).unwrap();
        assert!(read_trait_aligned(p.to_str().unwrap(),&names(&["a"])).is_err());
    }
    std::fs::remove_dir_all(base).unwrap();
}
