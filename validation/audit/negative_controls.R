args<-commandArgs(TRUE); root<-normalizePath(args[1]); base<-file.path(root,'validation/audit/artifacts/baseline'); out<-file.path(root,'validation/audit/artifacts/negative_controls');dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(root,'scripts/compare.R'))
results<-list()
for(test in c('unchanged','missing_file','numeric_change','unexpected_NA','unexpected_NaN','unexpected_Inf','unexpected_minusInf','duplicate_ID','unexpected_ID','reordered_ids','shuffled_labels','missing_colors','sample_ID_change','numeric_metrics','partition_renaming')) {
 dest<-file.path(out,test);dir.create(dest,showWarnings=FALSE);file.copy(list.files(file.path(base,'rust_out'),full.names=TRUE),dest,overwrite=TRUE)
 path<-file.path(dest,'cor.tsv');d<-read.delim(path,check.names=FALSE)
 if(test=='missing_file') unlink(path)
 if(test=='numeric_change') {d[1,2]<-d[1,2]+.1;write.table(d,path,sep='\t',quote=FALSE,row.names=FALSE)}
 if(test %in% c('unexpected_NA','unexpected_NaN','unexpected_Inf','unexpected_minusInf')) {d[1,2]<-switch(test,unexpected_NA=NA_real_,unexpected_NaN=NaN,unexpected_Inf=Inf,unexpected_minusInf=-Inf);write.table(d,path,sep='\t',quote=FALSE,row.names=FALSE)}
 if(test=='duplicate_ID') {d[2,1]<-d[1,1];write.table(d,path,sep='\t',quote=FALSE,row.names=FALSE)}
 if(test=='unexpected_ID') {d[1,1]<-'unknown_gene';write.table(d,path,sep='\t',quote=FALSE,row.names=FALSE)}
 if(test=='reordered_ids') for(f in c('adjacency.tsv','tom.tsv','eigengenes.tsv','kme.tsv','module_trait_cor.tsv','module_trait_p.tsv')) {
  x<-read.delim(file.path(dest,f),check.names=FALSE);x<-x[nrow(x):1,c(1,ncol(x):2)];write.table(x,file.path(dest,f),sep='\t',quote=FALSE,row.names=FALSE)
 }
 if(test=='shuffled_labels') {p<-file.path(dest,'module_colors.tsv');x<-readLines(p);writeLines(rev(x),p)}
 if(test=='missing_colors') unlink(file.path(dest,'module_colors.tsv'))
 if(test=='sample_ID_change') {p<-file.path(dest,'eigengenes.tsv');x<-read.delim(p,check.names=FALSE);x[1,1]<-'unknown_sample';write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE)}
 if(test=='numeric_metrics') {
  x<-matrix(c(0,2),2,dimnames=list(c('a','b'),'v'));y<-matrix(c(0,1),2,dimnames=dimnames(x));z<-error_metrics(x,y,.1)
  good<-z$max_abs==1 && z$mean_abs==.5 && abs(z$rmse-sqrt(.5))<1e-15 && z$above_tolerance==1 && z$total==2
  results[[test]]<-data.frame(test,expected_exit=0,observed_exit=0,status=if(good) 'PASS' else 'FAIL');next
 }
 if(test=='partition_renaming') {z<-partition_metrics(c('a','a','b','b'),c('x','x','y','y'));results[[test]]<-data.frame(test,expected_exit=0,observed_exit=0,status=if(z$ARI==1 && z$exact_agreement==0) 'PASS' else 'FAIL');next}
 cmd<-c(shQuote(file.path(root,'scripts/compare.R')),shQuote(dest),shQuote(file.path(base,'r_out')),'1e-9',shQuote(file.path(dest,'result')),shQuote(file.path(base,'testdata/expr.tsv')))
 code<-system2('Rscript',cmd,stdout=file.path(dest,'stdout.log'),stderr=file.path(dest,'stderr.log'))
 expected<-if(test %in% c('unchanged','reordered_ids')) 0 else 1
 results[[test]]<-data.frame(test,expected_exit=expected,observed_exit=code,status=if(code==expected) 'PASS' else 'FAIL')
}
r<-do.call(rbind,results);write.csv(r,file.path(out,'summary.csv'),row.names=FALSE);print(r);if(any(r$status!='PASS'))quit(status=1)
