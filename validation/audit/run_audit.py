import subprocess,json,time,pathlib,os,platform,hashlib
ROOT=pathlib.Path(__file__).resolve().parents[2]; A=ROOT/'validation/audit'; events=[]
def run(name,cmd,cwd=ROOT,env=None):
 t=time.time(); p=subprocess.run(cmd,cwd=cwd,shell=True,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=env)
 (A/'artifacts'/f'{name}.log').write_text(p.stdout)
 e=dict(name=name,command=cmd,cwd=str(cwd),exit_code=p.returncode,seconds=time.time()-t,log=f'artifacts/{name}.log');events.append(e)
 with (A/'test_execution_log.txt').open('a') as f:f.write(json.dumps(e)+'\n'+p.stdout+'\n')
 (A/'artifacts'/f'{name}.execution.json').write_text(json.dumps(e,indent=2));print(name,p.returncode,round(e['seconds'],2),flush=True)
 return p
if __name__=='__main__':
 run('environment',"rustc -Vv; cargo -V; Rscript -e 'sessionInfo(); print(extSoftVersion()); print(parallel::detectCores()); print(packageVersion(\"WGCNA\"))'; uname -a; sysctl -n machdep.cpu.brand_string hw.memsize hw.logicalcpu; git rev-parse HEAD; git ls-remote origin HEAD refs/heads/master")
 (A/'environment.json').write_text(json.dumps(dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),platform=platform.platform(),thread_environment={k:os.environ.get(k) for k in ['RAYON_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']},details='artifacts/environment.log'),indent=2))
 run('cargo_baseline','cargo test --workspace')
 b=A/'artifacts/baseline'; b.mkdir(exist_ok=True)
 for n in ['Cargo.toml','Cargo.lock','src','target']:
  (b/n).symlink_to(ROOT/n,target_is_directory=(ROOT/n).is_dir()) if not (b/n).exists() else None
 import shutil
 shutil.copytree(ROOT/'scripts',b/'scripts',dirs_exist_ok=True)
 run('parity_baseline',f'bash {b}/scripts/run_parity.sh 6 1e-9')
 run('build_real','cargo build --release --example real_bench')
 import sys
 sys.exit(1 if any(e['exit_code'] != 0 for e in events) else 0)
