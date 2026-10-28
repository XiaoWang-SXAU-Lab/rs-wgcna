from run_audit import *
import itertools,shlex,csv
import sys
if '--summarize-existing' in sys.argv:
 with (A/'parameter_coverage.csv').open() as f: rows=list(csv.DictReader(f))
 for row in rows:
  d=A/'artifacts'/f"matrix_{int(row['case']):03}"/'comparison'
  required=[d/'parity_details.csv',d/'extra_details.csv']
  good=all(p.exists() for p in required)
  for p in required:
   if p.exists():
    with p.open() as f: good=good and all(z['status']=='PASS' for z in csv.DictReader(f))
  row['status']='PASS' if good else 'FAIL'
 with (A/'parameter_coverage.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 print({s:sum(r['status']==s for r in rows) for s in ['PASS','FAIL']})
 sys.exit(1 if any(r['status']!='PASS' for r in rows) else 0)
configs=[]
for nt,tt,td in itertools.product(['unsigned','signed','signed hybrid'],['unsigned','signed','signed Nowick','unsigned 2','signed 2','signed Nowick 2'],['min','mean']):
 configs.append((nt,tt,td,6,2,.15,'standard'))
for p in [1,12,2.5]: configs.append(('unsigned','signed','min',p,2,.15,'standard'))
for ds,m in itertools.product([0,1,3,4],[0,.25,.5]): configs.append(('unsigned','signed','min',6,ds,m,'challenging'))
configs.append(('unsigned','none','min',6,2,.15,'standard'))
# Deterministic perturbation uses the existing synthetic fixture; no biological data changed.
run('generate_challenging',"Rscript -e 'x<-read.delim(\"validation/audit/artifacts/baseline/testdata/expr.tsv\",check.names=FALSE); set.seed(20261028); x[2:10,-1]<- -as.matrix(x[2:10,-1]); x[140:150,-1]<-as.matrix(x[140:150,-1])*.2+matrix(rnorm(11*40),11); write.table(x,\"validation/audit/artifacts/challenging.tsv\",sep=\"\\t\",quote=FALSE,row.names=FALSE)'")
rows=[]
for i,(nt,tt,td,p,ds,m,fixture) in enumerate(configs):
 dest=A/'artifacts'/f'matrix_{i:03}'
 if dest.exists(): raise RuntimeError(f'refusing overwrite: {dest}')
 dest.mkdir()
 base=A/'artifacts/baseline/testdata';expr=base/'expr.tsv' if fixture=='standard' else A/'artifacts/challenging.tsv'
 params=' '.join(map(shlex.quote,map(str,[p,base/'colors.txt',base/'trait.tsv',nt,tt,td,ds,m])))
 # Existing command argument order: expr out power colors trait nt tt td deep merge.
 rustargs=' '.join(map(shlex.quote,map(str,[expr,p,dest/'rust',base/'colors.txt',base/'trait.tsv',nt,tt,td,ds,m])))
 rargs=' '.join(map(shlex.quote,map(str,[expr,dest/'r',p,base/'colors.txt',base/'trait.tsv',nt,tt,td,ds,m])))
 pr=run(f'matrix_{i:03}_rust',f'target/release/examples/audit_matrix {rustargs}')
 rr=run(f'matrix_{i:03}_r',f'Rscript scripts/r_reference.R {rargs}')
 status='FAIL'
 if pr.returncode==0 and rr.returncode==0:
  pc=run(f'matrix_{i:03}_compare',f'Rscript scripts/compare.R {dest}/rust {dest}/r 1e-9 {dest}/comparison {expr}')
  pe=run(f'matrix_{i:03}_extra',f'Rscript validation/audit/compare_extra.R {dest}')
  status='PASS' if pc.returncode==0 and pe.returncode==0 else 'FAIL'
 rows.append(dict(case=i,network=nt,TOM=tt,denominator=td,power=p,deepSplit=ds,mergeCutHeight=m,fixture=fixture,status=status))
 with (A/'parameter_coverage.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

# A driver must not turn recorded required command failures into success.
import sys
sys.exit(1 if any(e["exit_code"] != 0 for e in events) else 0)
