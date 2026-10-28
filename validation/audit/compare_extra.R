source('scripts/compare.R');d<-commandArgs(TRUE)[1];rows<-list()
for(f in c('tree.tsv','dynamic.tsv','connectivity.tsv')) {
 x<-read_checked(file.path(d,'rust',f));y<-read_checked(file.path(d,'r',f));z<-error_metrics(x,y,1e-9)
 rows[[f]]<-cbind(metric=f,status=if(z$above_tolerance==0)'PASS' else 'FAIL',z)
 if(f=='dynamic.tsv') {p<-partition_metrics(as.vector(x),as.vector(y));write.csv(data.frame(exact_agreement=p$exact_agreement,ARI=p$ARI,disagreeing=p$disagreeing),file.path(d,'comparison/dynamic_summary.csv'),row.names=FALSE)}
}
r<-do.call(rbind,rows);write.csv(r,file.path(d,'comparison/extra_details.csv'),row.names=FALSE);print(r);if(any(r$status!='PASS'))quit(status=1)
