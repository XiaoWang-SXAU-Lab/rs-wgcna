from pathlib import Path
root=Path(__file__).resolve().parents[2]
s=(root/'examples/real_bench.rs').read_text()
s=s.replace('let _tom = tom_from_cor(', 'let tom = tom_from_cor(')
s=s.replace('    // 5. full module detection', '''    let probe: Vec<usize> = (0..g).step_by((g/64).max(1)).collect();
    let probe_genes: Vec<String> = probe.iter().map(|&i|filtered_genes[i].clone()).collect();
    for (name,mat) in [("cor_probe",&cor),("tom_probe",&tom)] {
        let sample = Mat::from_fn(probe.len(),g,|i,j|mat[(probe[i],j)]);
        write_mat(&format!("{outdir}/{name}.tsv"),&sample,&probe_genes,&filtered_genes);
    }
    let diss = Mat::from_fn(g,g,|i,j|1.0-tom[(i,j)]);
    let tree = rs_wgcna::hclust::hclust_average(&diss);
    let labels = rs_wgcna::cutree::cutree_hybrid(&tree.merge,&tree.height,&diss,30,2,0.995,true,true);
    let mut f = std::fs::File::create(format!("{outdir}/tree.tsv")).unwrap();
    writeln!(f,"ID\\tleft\\tright\\theight").unwrap();
    for (i,((a,b),h)) in tree.merge.iter().zip(tree.height.iter()).enumerate() {writeln!(f,"{}\\t{}\\t{}\\t{:.17}",i+1,a,b,h).unwrap();}
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
    // 5. full module detection''')
s=s.replace('    let traits = trait_data.1;', '''    assert_eq!(trait_data.0,exp.samples,"trait sample order differs");
    let traits = trait_data.1;''')
s=s.replace('    // write outputs', '''    write_mat(&format!("{outdir}/eigengenes.tsv"),&mes,&exp.samples,&mods);
    let km=rs_wgcna::kme::kme(&filtered,&mes);
    write_mat(&format!("{outdir}/kme.tsv"),&km,&filtered_genes,&mods);
    let intra=Mat::from_fn(g,1,|i,_|(0..g).filter(|&j|i!=j && colors[i]==colors[j]).map(|j|cor[(i,j)].abs().powi(6)).sum::<f64>());
    write_mat(&format!("{outdir}/intramodular.tsv"),&intra,&filtered_genes,&["within".into()]);
    // write outputs''')
s=s.replace('{:.6}', '{:.17}')
(root/'examples/audit_real.rs').write_text('// TEST-ONLY diagnostic extension of examples/real_bench.rs.\n'+s)
s=(root/'scripts/real_bench.R').read_text()
s=s.replace('datExpr <- t(expr)', '''stopifnot(!anyDuplicated(rownames(expr)), !anyDuplicated(colnames(expr)),all(is.finite(as.matrix(expr))))
datExpr <- t(expr)''')
s=s.replace('t <- proc.time()\nnet <- blockwiseModules', '''write_mat <- function(m,file) write.table(data.frame(ID=rownames(m),m,check.names=FALSE),file,sep="\\t",row.names=FALSE,quote=FALSE)
probe<-seq(1,ncol(datExpr),by=max(1,ncol(datExpr)%/%64));dimnames(tom)<-list(genes,genes)
write_mat(cor_mat[probe,,drop=FALSE],file.path(outdir,"cor_probe.tsv"));write_mat(tom[probe,,drop=FALSE],file.path(outdir,"tom_probe.tsv"))
tree<-hclust(as.dist(1-tom),method="average")
write.table(data.frame(ID=seq_along(tree$height),left=tree$merge[,1],right=tree$merge[,2],height=tree$height),file.path(outdir,"tree.tsv"),sep="\\t",quote=FALSE,row.names=FALSE)
dynamic<-cutreeDynamic(tree,distM=1-tom,deepSplit=2,minClusterSize=30,cutHeight=.995,pamRespectsDendro=TRUE,verbose=0)
write_mat(matrix(dynamic,ncol=1,dimnames=list(genes,"label")),file.path(outdir,"dynamic.tsv"))
sft<-pickSoftThreshold(datExpr,networkType="unsigned",powerVector=c(1:10,seq(12,20,2)),verbose=0)$fitIndices
write.table(sft,file.path(outdir,"sft.tsv"),sep="\\t",quote=FALSE,row.names=FALSE)
t <- proc.time()
net <- blockwiseModules''')
s=s.replace('trait_df <- read.table(traitpath, header = TRUE, row.names = 1, check.names = FALSE)', '''trait_df <- read.table(traitpath, header = TRUE, row.names = 1, check.names = FALSE)
stopifnot(!anyDuplicated(rownames(trait_df)),setequal(rownames(datExpr),rownames(trait_df)))''')
s=s.replace('# outputs (match Rust format)', '''# Export actual blockwiseModules premerge and final states. unmergedColors is
# after cleanup/reassignment, NOT the raw dynamic tree cut.
write.table(data.frame(ID=genes,color=net$unmergedColors),file.path(outdir,"premerge.tsv"),sep="\\t",quote=FALSE,row.names=FALSE)
colnames(MEs)<-substring(colnames(MEs),3)
write_mat(MEs,file.path(outdir,"eigengenes.tsv"))
km<-as.matrix(signedKME(datExpr,MEs));colnames(km)<-colnames(MEs)
write_mat(km,file.path(outdir,"kme.tsv"))
adj<-abs(cor_mat)^6;diag(adj)<-0
within<-vapply(seq_along(genes),function(i)sum(adj[i,net$colors==net$colors[i]]),numeric(1))
write_mat(matrix(within,ncol=1,dimnames=list(genes,"within")),file.path(outdir,"intramodular.tsv"));rm(adj)
# outputs (match Rust format)''')
s=s.replace('net <- blockwiseModules(', 'source("validation/audit/instrument_reference.R"); reference_diagnostic <- instrument_blockwise()\nnet <- reference_diagnostic(')
(root/'scripts/audit_real_reference.R').write_text('# TEST-ONLY diagnostic extension of scripts/real_bench.R.\n'+s)
