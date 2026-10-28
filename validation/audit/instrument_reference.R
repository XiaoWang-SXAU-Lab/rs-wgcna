# Diagnostic clone: record values after assignment, with expressions and their
# order unchanged. Original WGCNA namespace/function remains untouched.
instrument_blockwise <- function() {
 f <- WGCNA::blockwiseModules
 inject <- function(x) {
  if(!is.call(x)) return(x)
  if(identical(x[[1]],as.name('{'))) {
   out<-list(as.name('{'))
   for(i in seq.int(2,length(x))) {
    e<-inject(x[[i]]);out[[length(out)+1]]<-e
    if(is.call(e) && as.character(e[[1]])[1] %in% c('=','<-') && identical(e[[2]],as.name('goodLabels'))) {
      out[[length(out)+1]]<-quote(write.table(data.frame(ID=colnames(datExpr)[gsg$goodGenes],label=goodLabels),file.path(outdir,'cleanup.tsv'),sep='\t',quote=FALSE,row.names=FALSE))
    }
   }
   return(as.call(out))
  }
  for(i in seq_along(x)[-1]) if(is.call(x[[i]])) x[[i]]<-inject(x[[i]])
  x
 }
 body(f)<-inject(body(f)); e<-new.env(parent=environment(WGCNA::blockwiseModules))
 e$outdir<-outdir
 e$cutreeDynamic<-function(...) {
  z<-dynamicTreeCut::cutreeDynamic(...)
  write.table(data.frame(ID=genes,label=z),file.path(outdir,'dynamic_actual.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  z
 }
 environment(e$cutreeDynamic)<-environment()
 environment(f)<-e;f
}
