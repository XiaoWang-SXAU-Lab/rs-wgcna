from run_audit import *
b=A/'artifacts/diagnostic_smoke'; b.mkdir(exist_ok=True)
# Use original nonnegative biological fixture, top 150, unchanged parameters.
for impl,exe in [('rust','target/release/examples/audit_real'),('r','Rscript scripts/audit_real_reference.R')]:
 run(f'diagnostic_smoke_{impl}',f'{exe} Demo/WGCNA_fpkm_data.txt Demo/WGCNA_data_trait.txt {b}/{impl} 150')

# A driver must not turn recorded required command failures into success.
import sys
sys.exit(1 if any(e["exit_code"] != 0 for e in events) else 0)
