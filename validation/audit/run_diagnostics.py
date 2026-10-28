from run_audit import *
# Diagnostic runs add probes/intermediate exports: NOT performance measurements.
for n in [8000,12000]:
 for impl in ['rust','r']:
  dest=A/'artifacts'/f'diagnostic_{n}_{impl}'
  if dest.exists(): raise RuntimeError(f'refusing overwrite: {dest}')
  exe='target/release/examples/audit_real' if impl=='rust' else 'R_HOME=/Library/Frameworks/R.framework/Resources DYLD_LIBRARY_PATH='+str(A/'artifacts/accelerate_blas')+' VECLIB_MAXIMUM_THREADS=1 /Library/Frameworks/R.framework/Resources/bin/exec/R --no-save --no-restore --slave -f scripts/audit_real_reference.R --args'
  p=run(f'diagnostic_{n}_{impl}',f'{exe} Demo/WGCNA_fpkm_data.txt Demo/WGCNA_data_trait.txt {dest} {n}')

# A driver must not turn recorded required command failures into success.
import sys
sys.exit(1 if any(e["exit_code"] != 0 for e in events) else 0)
