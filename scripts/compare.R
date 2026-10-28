#!/usr/bin/env Rscript
# Strict parity comparison. No implicit optional files. Absolute tolerance decides
# PASS; relative error is diagnostic using max(abs(reference), 1e-12).
read_checked <- function(path) {
  if (!file.exists(path)) stop('required file missing: ', path)
  x <- read.delim(path, check.names=FALSE, row.names=NULL, stringsAsFactors=FALSE)
  if (ncol(x)<2 || nrow(x)<1) stop('empty/invalid table: ',path)
  ids <- as.character(x[[1]]); cols <- names(x)[-1]
  valid_ids <- function(z) length(z)>0 && !anyNA(z) && all(nzchar(z)) && !anyDuplicated(z)
  if (!valid_ids(ids) || !valid_ids(cols)) stop('missing/duplicate identifiers: ',path)
  m <- as.matrix(x[-1]); storage.mode(m)<-'double'; rownames(m)<-ids
  if (any(!is.finite(m))) stop('unexpected NA/NaN/Inf: ',path)
  m
}
align_checked <- function(x,y) {
  if (!identical(dim(x),dim(y))) stop('dimension mismatch')
  if (!setequal(rownames(x),rownames(y)) || !setequal(colnames(x),colnames(y))) stop('identifier set mismatch')
  x[rownames(y),colnames(y),drop=FALSE]
}
error_metrics <- function(x,y,tol) {
  x<-align_checked(x,y); d<-abs(x-y)
  data.frame(total=length(d),max_abs=max(d),mean_abs=mean(d),rmse=sqrt(mean(d*d)),
             max_relative=max(d/pmax(abs(y),1e-12)),above_tolerance=sum(d>tol),
             proportion_above=mean(d>tol),tolerance=tol)
}
partition_metrics <- function(x,y) {
  tab<-table(x,y); c2<-function(z) z*(z-1)/2; n<-length(x)
  a<-sum(c2(tab)); b<-sum(c2(rowSums(tab))); c<-sum(c2(colSums(tab)))
  expected<-if(n>1) b*c/c2(n) else 0; den<-(b+c)/2-expected
  ari<-if(den==0) if(a==(b+c)/2) 1 else 0 else (a-expected)/den
  list(total=n,exact_agreement=mean(x==y),disagreeing=sum(x!=y),ARI=ari,contingency=tab)
}
compare_main <- function(args) {
 if(length(args)<2) stop('usage: compare.R <rust_dir> <r_dir> [tol=1e-9] [results_dir] [expression.tsv]')
 tol<-if(length(args)>=3) as.numeric(args[3]) else 1e-9
 if(length(tol)!=1 || !is.finite(tol) || tol<0) stop('invalid tolerance')
 out<-if(length(args)>=4) args[4] else file.path(args[1],'comparison')
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 files<-c('cor.tsv','adjacency.tsv','tom.tsv','sft.tsv','eigengenes.tsv','kme.tsv','module_trait_cor.tsv','module_trait_p.tsv')
 rows<-list(); ok<-TRUE
 for(f in files) {
  row<-tryCatch({
   x<-read_checked(file.path(args[1],f)); y<-read_checked(file.path(args[2],f))
   if(f=='sft.tsv') { # Explicit WGCNA fitIndices naming aliases only.
    aliases<-c('SFT.R.sq'='R.sq','mean.k.'='mean.k','median.k.'='median.k','max.k.'='max.k')
    for(a in names(aliases)) {colnames(x)[colnames(x)==a]<-aliases[[a]];colnames(y)[colnames(y)==a]<-aliases[[a]]}
   }
   if(length(args)>=5) {
    ex<-read_checked(args[5]); expected_rows<-switch(f,'cor.tsv'=rownames(ex),'adjacency.tsv'=rownames(ex),'tom.tsv'=rownames(ex),'kme.tsv'=rownames(ex),'eigengenes.tsv'=colnames(ex),NULL)
    if(!is.null(expected_rows) && (!setequal(rownames(x),expected_rows)||!setequal(rownames(y),expected_rows))) stop('IDs differ from input expression')
    if(f %in% c('cor.tsv','adjacency.tsv','tom.tsv') && (!setequal(colnames(x),rownames(ex))||!setequal(colnames(y),rownames(ex)))) stop('column gene IDs differ from input expression')
   }
   z<-error_metrics(x,y,tol); cbind(metric=f,status=if(z$above_tolerance==0) 'PASS' else 'FAIL',reason='',z)
  },error=function(e) data.frame(metric=f,status='FAIL',reason=conditionMessage(e),total=NA,max_abs=NA,mean_abs=NA,rmse=NA,max_relative=NA,above_tolerance=NA,proportion_above=NA,tolerance=tol))
  rows[[f]]<-row; if(row$status!='PASS') ok<-FALSE
 }
 # Unlabelled color vectors require a proven gene-order contract from each
 # side's labelled correlation table; the generating harness uses expr order.
 tryCatch({
  gx<-rownames(read_checked(file.path(args[1],'cor.tsv'))); gy<-rownames(read_checked(file.path(args[2],'cor.tsv')))
  x<-readLines(file.path(args[1],'module_colors.tsv')); y<-readLines(file.path(args[2],'module_colors.tsv'))
  if(length(x)!=length(gx)||length(y)!=length(gy)||any(!nzchar(c(x,y)))) stop('invalid colors length/labels')
  if(!setequal(gx,gy)) stop('color gene IDs differ')
  if(length(args)>=5) {
   ex<-read_checked(args[5]); if(!identical(gx,rownames(ex))||!identical(gy,rownames(ex))) stop('unlabelled color ordering contract violated')
   meX<-read_checked(file.path(args[1],'eigengenes.tsv')); meY<-read_checked(file.path(args[2],'eigengenes.tsv'))
   if(!setequal(rownames(meX),colnames(ex))||!setequal(rownames(meY),colnames(ex))) stop('sample IDs differ from expression')
  }
  x<-x[match(gy,gx)]; z<-partition_metrics(x,y)
  write.csv(data.frame(total=z$total,exact_agreement=z$exact_agreement,disagreeing=z$disagreeing,ARI=z$ARI),file.path(out,'module_summary.csv'),row.names=FALSE)
  write.table(as.data.frame(z$contingency),file.path(out,'module_contingency.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(data.frame(gene=gy[x!=y],rust=x[x!=y],R=y[x!=y]),file.path(out,'module_disagreements.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  rows[['colors']]<-data.frame(metric='module_colors.tsv',status=if(z$disagreeing==0) 'PASS' else 'FAIL',reason=paste('exact disagreements',z$disagreeing,'ARI',z$ARI),total=z$total,max_abs=NA,mean_abs=NA,rmse=NA,max_relative=NA,above_tolerance=NA,proportion_above=NA,tolerance=tol)
  if(z$disagreeing>0) ok<-FALSE
 },error=function(e) {ok<<-FALSE; rows[['colors']]<<-data.frame(metric='module_colors.tsv',status='FAIL',reason=conditionMessage(e),total=NA,max_abs=NA,mean_abs=NA,rmse=NA,max_relative=NA,above_tolerance=NA,proportion_above=NA,tolerance=tol)})
 result<-do.call(rbind,rows);write.csv(result,file.path(out,'parity_details.csv'),row.names=FALSE);print(result,row.names=FALSE)
 cat('RESULT:',if(ok) 'PASS' else 'FAIL','\n'); if(!ok) quit(status=1)
}
if(sys.nframe()==0) compare_main(commandArgs(trailingOnly=TRUE))
