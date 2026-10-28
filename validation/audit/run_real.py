from run_audit import *
# Preserve each repeat in a new directory; refuse overwriting generated evidence.
for n in [8000,12000]:
 for rep in range(1, (3 if n==8000 else 1)+1):
  for impl in ['rust','r']:
   dest=A/'artifacts'/f'real_{n}_{impl}_{rep}'
   if dest.exists(): raise RuntimeError(f'refusing overwrite: {dest}')
   exe='target/release/examples/real_bench' if impl=='rust' else 'Rscript scripts/real_bench.R'
   run(f'real_{n}_{impl}_{rep}',f'/usr/bin/time -l {exe} Demo/WGCNA_fpkm_data.txt Demo/WGCNA_data_trait.txt {dest} {n}')

# A driver must not turn recorded required command failures into success.
import sys
sys.exit(1 if any(e["exit_code"] != 0 for e in events) else 0)
