suppressMessages(library(WGCNA));source('scripts/compare.R');d<-'validation/fixes/artifacts/stats_grid';dir.create(d,showWarnings=FALSE)
x<-read.delim(file.path(d,'rust.tsv'));y<-x;y$student<-mapply(corPvalueStudent,x$r,x$n);y$fisher<-mapply(corPvalueFisher,x$r,x$n)
write.table(y,file.path(d,'r.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
rows<-list();for(metric in c('student','fisher')) {
 a<-matrix(x[[metric]],ncol=1,dimnames=list(as.character(x$ID),metric));b<-matrix(y[[metric]],ncol=1,dimnames=dimnames(a));z<-error_metrics(a,b,1e-9)
 rows[[metric]]<-cbind(metric,status=if(z$above_tolerance==0)'PASS' else 'FAIL',z)
}
r<-do.call(rbind,rows);write.csv(r,file.path(d,'parity_details.csv'),row.names=FALSE);print(r)
z<-cbind(x,R_student=y$student,R_fisher=y$fisher,student_abs_error=abs(x$student-y$student),fisher_abs_error=abs(x$fisher-y$fisher));write.csv(z,file.path(d,'all_cases.csv'),row.names=FALSE)
print(z[order(z$student_abs_error,decreasing=TRUE)[1:6],]);if(any(r$status!='PASS'))quit(status=1)
