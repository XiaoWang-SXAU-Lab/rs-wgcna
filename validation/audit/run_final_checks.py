from run_audit import *
run('cargo_final','cargo test --workspace')
run('controls_final','Rscript validation/audit/negative_controls.R .')
run('parity_hardened','bash scripts/run_parity.sh 6 1e-9 validation/audit/artifacts/hardened')
# Trait ID reordering is scientifically legal. Existing Rust harness discards
# sample IDs; R aligns them. Preserve this additional minimal failing evidence.
base=A/'artifacts/baseline/testdata';d=A/'artifacts/trait_reordered';d.mkdir(exist_ok=True)
run('reorder_trait',f'''Rscript -e 'x<-read.delim("{base}/trait.tsv",check.names=FALSE);write.table(x[nrow(x):1,],"{d}/trait.tsv",sep="\\t",quote=FALSE,row.names=FALSE)' ''')
run('trait_reordered_rust',f'target/release/rs-wgcna {base}/expr.tsv 6 {d}/rust {base}/colors.txt {d}/trait.tsv')
run('trait_reordered_r',f'Rscript scripts/r_reference.R {base}/expr.tsv {d}/r 6 {base}/colors.txt {d}/trait.tsv')
run('trait_reordered_compare',f'Rscript scripts/compare.R {d}/rust {d}/r 1e-9 {d}/comparison {base}/expr.tsv')
run('input_alignment',"Rscript -e 'x<-read.delim(\"Demo/WGCNA_fpkm_data.txt\",check.names=FALSE,row.names=1);y<-read.delim(\"Demo/WGCNA_data_trait.txt\",check.names=FALSE,row.names=1);stopifnot(!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)),!anyDuplicated(rownames(y)),identical(colnames(x),rownames(y)),all(is.finite(as.matrix(x))),all(is.finite(as.matrix(y))));cat(nrow(x),ncol(x),\"gene/sample IDs unique and trait order identical\\n\")'")

# A driver must not turn recorded required command failures into success.
import sys
sys.exit(1 if any(e["exit_code"] != 0 for e in events) else 0)
