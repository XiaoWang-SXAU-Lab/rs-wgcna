//! Scientific parity regressions retained from the independent R audit.
use faer::Mat;
use rs_wgcna::network::{NetworkType,TomType,TomDenom};
fn fixture()->Mat<f64>{Mat::from_fn(2,2,|i,j|if i==j {1.0}else{0.25})}
#[test]
fn audit_fractional_adjacency_matches_reference() {
 let actual=rs_wgcna::adjacency::adjacency(&fixture(),2.5,NetworkType::Unsigned);
 assert!((actual[(0,1)]-0.03125).abs()<1e-12,"R |0.25|^2.5=0.03125; Rust={}",actual[(0,1)]);
}
#[test]
fn audit_fractional_tom_none_matches_reference() {
 let actual=rs_wgcna::tom::tom_from_cor(&fixture(),2.5,NetworkType::Unsigned,TomType::None,TomDenom::Min,1);
 assert!((actual[(0,1)]-0.03125).abs()<1e-12,"R TOMType none gives adjacency 0.03125; Rust={}",actual[(0,1)]);
}
#[test]
fn audit_near_zero_student_pvalue_matches_reference() {
 // WGCNA::corPvalueStudent(1e-8,5); two-sided p approaches one.
 let actual=rs_wgcna::stats::cor_pvalue(1e-8,5);
 let reference=0.9999999872676045_f64;
 assert!((actual-reference).abs()<1e-9,"R near-zero correlation p≈1; Rust={actual}");
}
#[test]
fn audit_signed_nowick2_unsigned_network_matches_reference() {
 // R TOMsimilarityFromExpr(cbind(as.numeric(1:4),-as.numeric(1:4)), networkType="unsigned",
 // TOMType="signed Nowick 2", power=6) returns off-diagonal +0.5.
 let cor=Mat::from_fn(2,2,|i,j|if i==j {1.0}else{-1.0});
 let actual=rs_wgcna::tom::tom_from_cor(&cor,6.0,NetworkType::Unsigned,TomType::SignedNowick2,TomDenom::Min,1);
 assert!((actual[(0,1)]-0.5).abs()<1e-9,"R Nowick2 TOM=+0.5; Rust={}",actual[(0,1)]);
}

#[test]
fn correlation_pvalue_endpoints_match_reference() {
    use rs_wgcna::stats::{cor_pvalue, cor_pvalue_fisher, normal_cdf};
    for r in [-1.0, 1.0] {
        assert_eq!(cor_pvalue(r, 56), 0.0);
        assert_eq!(cor_pvalue_fisher(r, 56), 0.0);
    }
    assert_eq!(cor_pvalue(0.0, 56), 1.0);
    assert_eq!(cor_pvalue_fisher(0.0, 56), 1.0);
    assert_eq!(normal_cdf(f64::NEG_INFINITY), 0.0);
    assert_eq!(normal_cdf(f64::INFINITY), 1.0);
}
