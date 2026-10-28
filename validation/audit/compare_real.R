source('scripts/compare.R')
args<-commandArgs(TRUE);xdir<-args[1];ydir<-args[2];out<-args[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_labels<-function(d){x<-read.delim(file.path(d,'module_colors.txt'),header=FALSE,stringsAsFactors=FALSE);names(x)<-c('gene','color');if(anyDuplicated(x$gene)||anyNA(x)||any(!nzchar(x$gene)))stop('invalid labels');x}
x<-read_labels(xdir);y<-read_labels(ydir);stopifnot(setequal(x$gene,y$gene));order_equal<-identical(x$gene,y$gene);x<-x[match(y$gene,x$gene),]
genes<-readLines(file.path(xdir,'filtered_genes.txt'));refgenes<-readLines(file.path(ydir,'filtered_genes.txt'))
stopifnot(!anyDuplicated(genes),!anyDuplicated(refgenes),identical(genes,x$gene),identical(refgenes,y$gene))
p<-partition_metrics(x$color,y$color)
write.csv(data.frame(total=p$total,exact_agreement=p$exact_agreement,disagreeing=p$disagreeing,ARI=p$ARI,gene_order_identical=order_equal,status=if(p$disagreeing==0)'PASS' else 'FAIL'),file.path(out,'module_summary.csv'),row.names=FALSE)
write.table(as.data.frame(p$contingency),file.path(out,'module_contingency.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
dis<-data.frame(gene=y$gene[x$color!=y$color],rust=x$color[x$color!=y$color],R=y$color[x$color!=y$color])
write.table(dis,file.path(out,'module_disagreements.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
rows<-list()
# Six-decimal historical/default benchmark files cannot validate machine precision.
# Strict 1e-9 remains the required threshold. Quantization-aware 5e-7 is an
# additional reporting bound only, never substituted into PASS decisions.
for(f in c('module_trait_cor.tsv','module_trait_p.tsv','cor_probe.tsv','tom_probe.tsv','sft.tsv','tree.tsv','dynamic.tsv','eigengenes.tsv','kme.tsv','intramodular.tsv')) {
 if(!file.exists(file.path(xdir,f))||!file.exists(file.path(ydir,f))) {rows[[f]]<-data.frame(metric=f,status=if(grepl('diagnostic',xdir)||f %in% c('module_trait_cor.tsv','module_trait_p.tsv'))'FAIL' else 'NOT EXECUTED',reason='required diagnostic file missing or not exported by original benchmark',total=NA,max_abs=NA,mean_abs=NA,rmse=NA,max_relative=NA,above_tolerance=NA,proportion_above=NA,tolerance=1e-9);next}
 rows[[f]]<-tryCatch({
  a<-read_checked(file.path(xdir,f));b<-read_checked(file.path(ydir,f));notes<-''
  if(f %in% c('module_trait_cor.tsv','module_trait_p.tsv')) {rownames(a)<-sub('^ME','',rownames(a));rownames(b)<-sub('^ME','',rownames(b))}
  if(f=='sft.tsv'){aliases<-c('SFT.R.sq'='R.sq','mean.k.'='mean.k','median.k.'='median.k','max.k.'='max.k');for(v in names(aliases))colnames(b)[colnames(b)==v]<-aliases[[v]]}
  if(f %in% c('module_trait_cor.tsv','module_trait_p.tsv','eigengenes.tsv','kme.tsv')) {
   # R's grey eigengene is expressly excluded by Rust real_bench.rs. Compare
   # the non-grey scientific estimands and report exclusion, not silent drops.
   if(f %in% c('eigengenes.tsv','kme.tsv')){b<-b[,colnames(b)!='grey',drop=FALSE]}else b<-b[rownames(b)!='grey',,drop=FALSE]
   notes<-'R grey excluded explicitly: Rust benchmark exports non-grey only; final-module differences may propagate'
  }
  z<-error_metrics(a,b,1e-9);cbind(metric=f,status=if(z$above_tolerance==0)'PASS' else 'FAIL',reason=notes,z)
 },error=function(e)data.frame(metric=f,status='FAIL',reason=conditionMessage(e),total=NA,max_abs=NA,mean_abs=NA,rmse=NA,max_relative=NA,above_tolerance=NA,proportion_above=NA,tolerance=1e-9))
}
write.csv(do.call(rbind,rows),file.path(out,'parity_details.csv'),row.names=FALSE)
if(file.exists(file.path(xdir,'dynamic.tsv')) && file.exists(file.path(ydir,'dynamic.tsv'))) {
 a<-read_checked(file.path(xdir,'dynamic.tsv'));b<-read_checked(file.path(ydir,'dynamic.tsv'));a<-align_checked(a,b)
 z<-partition_metrics(a[,1],b[,1]);write.csv(data.frame(total=z$total,exact_agreement=z$exact_agreement,disagreeing=z$disagreeing,ARI=z$ARI),file.path(out,'dynamic_summary.csv'),row.names=FALSE)
}
# Annotate every discordant gene using genuinely exported stages and kME.
if(nrow(dis)>0 && file.exists(file.path(xdir,'kme.tsv'))) {
 fetch<-function(dir,file,id){p<-file.path(dir,file);if(!file.exists(p))return(NA_character_);z<-read.delim(p,check.names=FALSE);as.character(z[match(id,z[[1]]),2])}
 a<-read_checked(file.path(xdir,'kme.tsv'));b<-read_checked(file.path(ydir,'kme.tsv'))
 for(i in seq_len(nrow(dis))){g<-dis$gene[i];dis$rust_dynamic[i]<-fetch(xdir,'dynamic.tsv',g);dis$R_dynamic[i]<-fetch(ydir,'dynamic.tsv',g);dis$rust_cleanup[i]<-fetch(xdir,'cleanup.tsv',g);dis$R_cleanup[i]<-fetch(ydir,'cleanup.tsv',g);dis$R_premerge[i]<-fetch(ydir,'premerge.tsv',g)
  for(m in union(colnames(a),colnames(b))) {dis[i,paste0('rust_kME_',m)]<-if(m %in% colnames(a))a[g,m] else NA;dis[i,paste0('R_kME_',m)]<-if(m %in% colnames(b))b[g,m] else NA}
  dis$rust_final_own_kME[i]<-a[g,dis$rust[i]];dis$R_final_own_kME[i]<-b[g,dis$R[i]]
  dis$rust_final_distance_minKMEtoStay[i]<-abs(a[g,dis$rust[i]])-.3;dis$R_final_distance_minKMEtoStay[i]<-abs(b[g,dis$R[i]])-.3
 }
 write.table(dis,file.path(out,'mismatch_investigation.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
}
print(read.csv(file.path(out,'module_summary.csv')));print(do.call(rbind,rows));if(p$disagreeing>0 || any(vapply(rows,function(x)x$status=='FAIL',logical(1))))quit(status=1)
