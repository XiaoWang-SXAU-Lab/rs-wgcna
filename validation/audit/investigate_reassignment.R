suppressMessages(library(WGCNA));source('scripts/compare.R')
base<-'validation/audit/artifacts';d<-file.path(base,'diagnostic_12000_r');u<-file.path(base,'diagnostic_12000_rust');genes<-readLines(file.path(d,'filtered_genes.txt'))
e<-read.delim('Demo/WGCNA_fpkm_data.txt',row.names=1,check.names=FALSE);x<-t(log10(as.matrix(e[genes,,drop=FALSE])+1));g<-'ENSG00000149541'
read_stage<-function(dir,f){z<-read.delim(file.path(dir,f),check.names=FALSE);stopifnot(setequal(z[[1]],genes));as.character(z[match(genes,z[[1]]),2])}
rows<-list()
for(impl in c('r','rust'))for(stage in c('dynamic','cleanup')) {
 dir<-if(impl=='r')d else u;labels<-read_stage(dir,paste0(stage,'.tsv'));ids<-sort(setdiff(unique(labels),'0'))
 me<-moduleEigengenes(x[,labels!='0',drop=FALSE],labels[labels!='0'])$eigengenes
 k<-abs(as.numeric(cor(x[,g],me)));names(k)<-substring(colnames(me),3);own<-labels[match(g,genes)];best<-names(which.max(k));ratio<-corPvalueFisher(k[best],nrow(x))/corPvalueFisher(k[own],nrow(x))
 rows[[paste(impl,stage)]]<-data.frame(impl,stage,gene=g,own_label=own,best_label=best,own_abs_kME=k[own],best_abs_kME=k[best],p_own=corPvalueFisher(k[own],nrow(x)),p_best=corPvalueFisher(k[best],nrow(x)),ratio=ratio,threshold=1e-6,distance=ratio-1e-6,reassign=ratio<1e-6)
}
r<-do.call(rbind,rows);write.table(r,file.path(base,'diagnostic_12000_comparison/reassignment_thresholds.tsv'),sep='\t',quote=FALSE,row.names=FALSE);print(r)
# Determine if actual blockwise cut equals independently reconstructed cut.
a<-read_stage(d,'dynamic_actual.tsv');b<-read_stage(d,'dynamic.tsv');print(table(actual=a,reconstructed=b));stopifnot(identical(a,b))
writeLines('PASS: actual R blockwise cut labels identical to standalone exported dynamic labels',file.path(base,'diagnostic_12000_comparison/actual_cut_check.txt'))
