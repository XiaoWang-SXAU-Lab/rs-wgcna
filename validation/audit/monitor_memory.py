import subprocess,time,os,json,pathlib
from run_audit import A
pid=34401;observations={}
while True:
 s=subprocess.run(['ps','-axo','pid=,rss=,stat=,command='],text=True,stdout=subprocess.PIPE).stdout
 live=False
 for line in s.splitlines():
  bits=line.strip().split(None,3)
  if len(bits)<4:continue
  p,rss,state,cmd=bits
  if int(p)==pid and not state.startswith('Z'):live=True
  if ('target/release/examples/audit_perf Demo/' in cmd or '--slave -f scripts/real_bench.R --args Demo/' in cmd) and not cmd.startswith('/bin/zsh'):
   k=p+' '+cmd;z=observations.setdefault(k,dict(pid=int(p),command=cmd,max_sampled_rss_bytes=0,samples=0));z['max_sampled_rss_bytes']=max(z['max_sampled_rss_bytes'],int(rss)*1024);z['samples']+=1
 (A/'artifacts/perfB_sampled_memory.json').write_text(json.dumps(dict(note='1 second ps RSS samples; lower-bound observations, not exact process peak; monitoring started after first Rust repeat',processes=list(observations.values())),indent=2))
 if not live:break
 time.sleep(1)
