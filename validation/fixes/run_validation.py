"""Validate fixes with the existing audit fixtures/comparators, preserving old evidence."""
import subprocess,pathlib,json,time,csv,shlex,os,hashlib,argparse
ROOT=pathlib.Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir',type=pathlib.Path,default=ROOT/'validation/fixes',help='fresh result directory; existing evidence is never overwritten')
options=parser.parse_args()
OUT=options.output_dir.resolve();ART=OUT/'artifacts';events=[]
def run(name,args,env=None):
    t=time.time();p=subprocess.run(args,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=env)
    (ART/f'{name}.log').write_text(p.stdout)
    record=dict(name=name,argv=[str(a) for a in args],cwd=str(ROOT),exit_code=p.returncode,seconds=time.time()-t,log=f'artifacts/{name}.log')
    events.append(record);(OUT/'test_manifest.json').write_text(json.dumps(events,indent=2));print(name,p.returncode,round(record['seconds'],2),flush=True)
    return p.returncode

def main():
    if (OUT/'test_manifest.json').exists():
        raise RuntimeError('Existing fix evidence: use fresh output paths with the recorded commands; refusing overwrite.')
    ART.mkdir(parents=True,exist_ok=True)
    assert run('cargo_tests',['cargo','test','--workspace'])==0
    assert run('build_release',['cargo','build','--release','--bins','--example','audit_matrix','--example','audit_real','--example','audit_stats','--example','real_bench'])==0
    base=ROOT/'validation/audit/artifacts/baseline/testdata'
    assert run('core_parity',['bash','scripts/run_parity.sh','6','1e-9',str(ART/'core')])==0
    # Old controls mutate preserved pre-fix outputs; this checks the validator,
    # not whether old numerical output was regenerated to make controls pass.
    controls=ROOT/'validation/audit/negative_controls.R'
    cs=controls.read_text().replace("file.path(root,'validation/audit/artifacts/negative_controls')",json.dumps(str(ART/'negative_controls')))
    (OUT/'negative_controls.R').write_text(cs)
    assert run('validator_controls',['Rscript',str(OUT/'negative_controls.R'),str(ROOT)])==0
    # The stats reference script only changes its output directory, not its formulas/tolerance.
    stats=ART/'stats_grid';stats.mkdir()
    proc=subprocess.run([str(ROOT/'target/release/examples/audit_stats')],cwd=ROOT,text=True,stdout=subprocess.PIPE,check=True)
    (stats/'rust.tsv').write_text(proc.stdout)
    s=(ROOT/'validation/audit/stats_reference.R').read_text().replace("d<-'validation/audit/artifacts/stats_grid'","d<-"+json.dumps(str(stats)))
    (OUT/'stats_reference.R').write_text(s)
    assert run('stats_reference',['Rscript',str(OUT/'stats_reference.R')])==0
    # Explicitly reverse trait rows using the same audited fixture and compare against R.
    trait=ROOT/'validation/audit/artifacts/trait_reordered/trait.tsv';dest=ART/'trait_reordered';dest.mkdir()
    assert run('trait_reordered_rust',[str(ROOT/'target/release/rs-wgcna'),str(base/'expr.tsv'),'6',str(dest/'rust'),str(base/'colors.txt'),str(trait)])==0
    assert run('trait_reordered_R',['Rscript','scripts/r_reference.R',str(base/'expr.tsv'),str(dest/'r'),'6',str(base/'colors.txt'),str(trait)])==0
    assert run('trait_reordered_compare',['Rscript','scripts/compare.R',str(dest/'rust'),str(dest/'r'),'1e-9',str(dest/'comparison'),str(base/'expr.tsv')])==0
    coverage=list(csv.DictReader((ROOT/'validation/audit/parameter_coverage.csv').open()));new=[]
    for row in coverage:
        d=ART/f"matrix_{int(row['case']):03}";d.mkdir()
        expr=base/'expr.tsv' if row['fixture']=='standard' else ROOT/'validation/audit/artifacts/challenging.tsv'
        config=[row[k] for k in ['network','TOM','denominator','deepSplit','mergeCutHeight']]
        common=[str(base/'colors.txt'),str(base/'trait.tsv')]+config
        c1=run(d.name+'_rust',[str(ROOT/'target/release/examples/audit_matrix'),str(expr),row['power'],str(d/'rust')]+common)
        c2=run(d.name+'_R',['Rscript','scripts/r_reference.R',str(expr),str(d/'r'),row['power']]+common)
        c3=run(d.name+'_compare',['Rscript','scripts/compare.R',str(d/'rust'),str(d/'r'),'1e-9',str(d/'comparison'),str(expr)]) if c1==0 and c2==0 else 1
        c4=run(d.name+'_extra',['Rscript','validation/audit/compare_extra.R',str(d)]) if c1==0 and c2==0 else 1
        new.append(dict(row,before_status=row['status'],status='PASS' if not any([c1,c2,c3,c4]) else 'FAIL'))
        with (OUT/'parameter_coverage.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=new[0]);w.writeheader();w.writerows(new)
    # Real-data diagnostics, not timing claims. Reference uses the already
    # verified process-local Accelerate BLAS. Production defaults stay unchanged.
    # Portable default: use the installed Rscript and its configured BLAS.
    # The archived audit used process-local Accelerate, recorded in its logs.
    env=dict(os.environ)
    for n in [8000,12000]:
        for impl in ['rust','r']:
            dest=ART/f'diagnostic_{n}_{impl}'
            exe=[str(ROOT/'target/release/examples/audit_real')] if impl=='rust' else ['Rscript','scripts/audit_real_reference.R']
            code=run(f'diagnostic_{n}_{impl}',exe+['Demo/WGCNA_fpkm_data.txt','Demo/WGCNA_data_trait.txt',str(dest),str(n)],env=env if impl=='r' else None)
            assert code==0
        run(f'diagnostic_{n}_compare',['Rscript','validation/audit/compare_real.R',str(ART/f'diagnostic_{n}_rust'),str(ART/f'diagnostic_{n}_r'),str(ART/f'diagnostic_{n}_comparison')])
    return 1 if any(e['exit_code'] for e in events) else 0
if __name__=='__main__':raise SystemExit(main())
