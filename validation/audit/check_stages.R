source('scripts/compare.R');rows<-list()
for(n in c(8000,12000)) {
 b<-'validation/audit/artifacts';u<-file.path(b,paste0('diagnostic_',n,'_rust'));r<-file.path(b,paste0('diagnostic_',n,'_r'))
 for(stage in c('dynamic','cleanup')){
  x<-read_checked(file.path(u,paste0(stage,'.tsv')));y<-read_checked(file.path(r,paste0(stage,'.tsv')));x<-align_checked(x,y);z<-partition_metrics(x[,1],y[,1])
  rows[[paste(n,stage)]]<-data.frame(genes=n,stage,total=z$total,exact_agreement=z$exact_agreement,ARI=z$ARI,disagreeing=z$disagreeing,status=if(z$disagreeing==0)'PASS' else 'FAIL')
 }
 x<-read_checked(file.path(r,'dynamic_actual.tsv'));y<-read_checked(file.path(r,'dynamic.tsv'));x<-align_checked(x,y);z<-partition_metrics(x[,1],y[,1]);rows[[paste(n,'actual_cut')]]<-data.frame(genes=n,stage='actual_R_cut_vs_reconstructed',total=z$total,exact_agreement=z$exact_agreement,ARI=z$ARI,disagreeing=z$disagreeing,status=if(z$disagreeing==0)'PASS' else 'FAIL')
}
z<-do.call(rbind,rows);write.csv(z,'validation/audit/stage_summary.csv',row.names=FALSE);print(z);if(any(z$status!='PASS'))quit(status=1)
