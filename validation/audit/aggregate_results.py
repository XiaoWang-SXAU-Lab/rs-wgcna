from run_audit import ROOT,A
import csv,json,pathlib,re,hashlib,statistics,subprocess,platform,os

def read(p):
 with open(p,newline='') as f:return list(csv.DictReader(f))
def write(p,rows,fields=None):
 if fields is None:fields=list(dict.fromkeys(k for r in rows for k in r))
 with open(p,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
# Keep every command attempt (including recovered failures), not merely latest.
events=[]
for line in (A/'test_execution_log.txt').read_text().splitlines():
 if line.startswith('{"name":'):
  try:z=json.loads(line)
  except ValueError:continue
  if 'command' in z and 'exit_code' in z:events.append(z)
(A/'test_manifest.json').write_text(json.dumps(events,indent=2,ensure_ascii=False))
summary=[];details=[];dis=[]
paths=[A/'artifacts/hardened_final/comparison']+list((A/'artifacts').glob('matrix_*/comparison'))+list((A/'artifacts').glob('*_comparison*'))+[A/'artifacts/trait_reordered/comparison']
for d in paths:
 if not d.is_dir():continue
 if (d/'parity_details.csv').exists():
  for z in read(d/'parity_details.csv'):details.append(dict(evidence=str((d/'parity_details.csv').relative_to(ROOT)),case=str(d.relative_to(A)),**z))
 if (d/'extra_details.csv').exists():
  for z in read(d/'extra_details.csv'):details.append(dict(evidence=str((d/'extra_details.csv').relative_to(ROOT)),case=str(d.relative_to(A)),**z))
 for name in ['module','dynamic']:
  f=d/f'{name}_summary.csv'
  if f.exists():
   for z in read(f):summary.append(dict(case=str(d.relative_to(A)),kind=name,evidence=str(f.relative_to(ROOT)),**z))
 f=d/'module_disagreements.tsv'
 if f.exists():
  with f.open() as h:
   for z in csv.DictReader(h,delimiter='\t'):dis.append(dict(case=str(d.relative_to(A)),evidence=str(f.relative_to(ROOT)),**z))
f=A/'artifacts/stats_grid/parity_details.csv'
if f.exists():
 for z in read(f):details.append(dict(evidence=str(f.relative_to(ROOT)),case='stats_grid',**z))
write(A/'parity_details.csv',details);write(A/'parity_summary.csv',summary)
with (A/'module_disagreements.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['case','evidence','gene','rust','R'],delimiter='\t');w.writeheader();w.writerows(dis)
# Performance: no cor/TOM double count in non-overlapping total. Preserve per
# repeat data and annotate default timings that overlapped audit work.
raw=[]
for folder in (A/'artifacts').glob('*_*_*_*'):
 m=re.fullmatch(r'(real|perfB)_(8000|12000)_(rust|r)_([123])',folder.name)
 if not m or not (folder/'timings.tsv').exists():continue
 with (folder/'timings.tsv').open() as f:ts=list(csv.DictReader(f,delimiter='\t'))
 vals={v['step']:float(v['seconds']) for v in ts}
 vals['nonoverlap_total']=sum(vals[k] for k in ['log10_normalize','mad_topN','module_detection_total','module_trait'])
 for step,sec in vals.items():raw.append(dict(scenario='A' if m[1]=='real' else 'B',genes=int(m[2]),implementation=m[3],repeat=int(m[4]),step=step,seconds=sec,evidence=str((folder/'timings.tsv').relative_to(ROOT))))
write(A/'benchmark_details.csv',raw)
groups={}
for r in raw:groups.setdefault((r['scenario'],r['genes'],r['implementation'],r['step']),[]).append(r['seconds'])
bench=[]
for (scenario,n,impl,step),v in sorted(groups.items()):
 bench.append(dict(scenario=scenario,genes=n,implementation=impl,step=step,repeats=len(v),median_seconds=statistics.median(v),min_seconds=min(v),max_seconds=max(v),sd_seconds=statistics.stdev(v) if len(v)>1 else '',threads=8 if scenario=='B' or impl=='rust' else 1,BLAS='faer' if impl=='rust' else ('reference libRblas' if scenario=='A' else 'Accelerate'),measurement_note='audit workloads overlapped some repeats; descriptive only' if scenario=='A' else 'sequential controlled-thread audit runs; other desktop activity uncontrolled'))
write(A/'benchmark_summary.csv',bench)
# Machine-level environment and file hashes.
env=json.loads((A/'environment.json').read_text());env.update(cpu='Apple M2',logical_cores=8,ram_bytes=17179869184,rustc='1.99.1 (58cf8f9c3 2026-10-13)',cargo='1.99.1',default_R=json.loads((A/'artifacts/R_environment_default.json').read_text()),scenario_B_R=json.loads((A/'artifacts/R_environment_accelerate8.json').read_text()),diagnostic_R=dict(BLAS='Accelerate',VECLIB_MAXIMUM_THREADS=1,WGCNA_threads=1),date_local='2026-10-28',timezone='Asia/Shanghai')
(A/'environment.json').write_text(json.dumps(env,indent=2,ensure_ascii=False))
hashes={}
for pat in ['src/**/*.rs','Demo/*','bench_out/**/*','scripts/*','examples/audit*.rs','tests/audit*.rs','report/benchmark.qmd','Cargo.toml','Cargo.lock']:
 for p in ROOT.glob(pat):
  if p.is_file():hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
(A/'file_hashes.json').write_text(json.dumps(hashes,indent=2))
checks=[]
def count(group,statuses):
 checks.append(dict(group=group,**{s:statuses.count(s) for s in ['PASS','FAIL','SKIPPED','NOT SUPPORTED','NOT EXECUTED']}))
count('original Rust unit baseline',['PASS']*12)
count('original synthetic baseline',['PASS']*9)
# Actual final Rust regression counts come from executable log, not intent.
f=A/'artifacts/cargo_verified_final.log'
if f.exists():
 statuses=[]
 for p,q in re.findall(r'test result: (?:ok|FAILED)\. (\d+) passed; (\d+) failed;',f.read_text()):statuses+=['PASS']*int(p)+['FAIL']*int(q)
 count('final cargo tests (original + new regressions)',statuses)
f=A/'parameter_coverage.csv'
if f.exists():count('parameter configurations',[x['status'] for x in read(f)])
f=A/'artifacts/negative_controls/summary.csv'
if f.exists():count('validator controls including expected failures',[x['status'] for x in read(f)])
count('comparison rows current artifacts (historical excluded)',[x['status'] for x in details if 'historical' not in x['case']])
count('statistical grid cases Student',[('PASS' if float(x['student_abs_error'])<=1e-9 else 'FAIL') for x in read(A/'artifacts/stats_grid/all_cases.csv')])
count('statistical grid cases Fisher',[('PASS' if float(x['fisher_abs_error'])<=1e-9 else 'FAIL') for x in read(A/'artifacts/stats_grid/all_cases.csv')])
count('additional default 12K repetitions omitted by runtime policy',['SKIPPED']*2)
f=A/'remote_final.json'
if f.exists():count('end-of-audit remote recheck',[json.loads(f.read_text())['status']])
write(A/'test_counts.csv',checks)
print('aggregated',len(events),'commands',len(details),'comparison rows',len(summary),'partition summaries')
