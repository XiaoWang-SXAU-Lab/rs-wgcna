from run_audit import *
# Wait for original default-environment measurements and diagnostics to finish.
# This permits independent comparison of quiet, controlled Scenario B.
import sys
for pid in map(int,sys.argv[1:]):
 while True:
  try: os.kill(pid,0)
  except ProcessLookupError: break
  time.sleep(3)
for n in [8000,12000]:
 for rep in [1,2,3]:
  for impl in ['rust','r']:
   dest=A/'artifacts'/f'perfB_{n}_{impl}_{rep}'
   if dest.exists():raise RuntimeError('refusing overwrite '+str(dest))
   if impl=='rust': exe='AUDIT_THREADS=8 RAYON_NUM_THREADS=8 target/release/examples/audit_perf'
   else:exe='R_HOME=/Library/Frameworks/R.framework/Resources DYLD_LIBRARY_PATH='+str(A/'artifacts/accelerate_blas')+' VECLIB_MAXIMUM_THREADS=8 ALLOW_WGCNA_THREADS=8 /Library/Frameworks/R.framework/Resources/bin/exec/R --no-save --no-restore --slave -f scripts/real_bench.R --args'
   run(f'perfB_{n}_{impl}_{rep}',f'{exe} Demo/WGCNA_fpkm_data.txt Demo/WGCNA_data_trait.txt {dest} {n}')

# A driver must not turn recorded required command failures into success.
import sys
sys.exit(1 if any(e["exit_code"] != 0 for e in events) else 0)
