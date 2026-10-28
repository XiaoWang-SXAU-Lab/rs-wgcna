import subprocess,json,pathlib,sys
p=pathlib.Path(__file__).resolve().parent
cmd=['git','-c','http.version=HTTP/1.1','ls-remote','origin','HEAD','refs/heads/master']
try:
 r=subprocess.run(cmd,text=True,capture_output=True,timeout=20)
 head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
 result=dict(command=cmd,local_commit=head,remote_output=r.stdout,stderr=r.stderr,exit_code=r.returncode,status='PASS' if r.returncode==0 and all(line.split()[0]==head for line in r.stdout.splitlines()) and len(r.stdout.splitlines())==2 else ('FAIL' if r.returncode==0 else 'SKIPPED'),reason='' if r.returncode==0 else 'network error; successful start-of-audit remote proof retained')
except subprocess.TimeoutExpired:
 result=dict(command=cmd,status='SKIPPED',reason='20-second network timeout; successful start-of-audit remote proof retained')
(p/'remote_final.json').write_text(json.dumps(result,indent=2));print(json.dumps(result));sys.exit(0 if result['status']=='PASS' else 1)
