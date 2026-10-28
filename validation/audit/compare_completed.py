from run_audit import *
for n in [8000,12000]:
 for prefix in ['real','perfB']:
  for rep in [1,2,3]:
   x=A/'artifacts'/f'{prefix}_{n}_rust_{rep}';y=A/'artifacts'/f'{prefix}_{n}_r_{rep}';dest=A/'artifacts'/f'{prefix}_{n}_comparison_{rep}'
   if (x/'timings.tsv').exists() and (y/'timings.tsv').exists() and not (dest/'parity_details.csv').exists():
    run(f'{prefix}_{n}_compare_{rep}',f'Rscript validation/audit/compare_real.R {x} {y} {dest}')
 for name,x,y in [('historical',ROOT/'bench_out'/('rust12k' if n==12000 else 'rust'),ROOT/'bench_out'/('r12k' if n==12000 else 'r')),('diagnostic',A/'artifacts'/f'diagnostic_{n}_rust',A/'artifacts'/f'diagnostic_{n}_r')]:
  if (A/'artifacts'/f'{name}_{n}_comparison'/'parity_details.csv').exists(): continue
  run(f'{name}_{n}_compare_verified',f'Rscript validation/audit/compare_real.R {x} {y} {A}/artifacts/{name}_{n}_comparison')

# A driver must not turn recorded required command failures into success.
import sys
sys.exit(1 if any(e["exit_code"] != 0 for e in events) else 0)
