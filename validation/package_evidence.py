"""List compact review evidence without deleting or changing local raw exports.

This builds a checksum manifest and prints artifact paths suitable for explicit
`git add -f`. No files are staged automatically; symlinks are never selected.
"""
from pathlib import Path
import hashlib,json,re
ROOT=Path(__file__).resolve().parents[1]

def selected(p):
    rel=p.relative_to(ROOT);name=p.name;parts=rel.parts
    if p.is_symlink() or '__pycache__' in parts:return False
    if 'artifacts' not in parts:
        return name!='FIX_VERIFICATION.generated.zh-CN.md' and p.suffix in {'.md','.csv','.tsv','.json','.txt','.py','.R','.patch'}
    sub=parts[parts.index('artifacts')+1:]
    if not sub:return False
    # Log records and command metadata are small, reviewable execution evidence.
    if len(sub)==1 and p.suffix in {'.log','.json','.R','.qmd'}:return True
    if any(re.fullmatch(r'(?:.*_)?comparison(?:_\d+)?',x) for x in sub[:-1]):
        return p.suffix in {'.csv','.tsv','.txt'}
    if name in {'summary.csv','stdout.log','stderr.log'}:return True
    # Preserve original control inputs and fresh default-core outputs in full.
    if sub[0]=='baseline' and len(sub)>1 and sub[1] in {'scripts','testdata','r_out','rust_out'}:return True
    if sub[0]=='core' and len(sub)>1 and sub[1] in {'testdata','r_out','rust_out','comparison'}:return True
    if sub[0]=='stats_grid':return p.suffix in {'.csv','.tsv'}
    if tuple(sub)==('challenging.tsv',) or tuple(sub)==('trait_reordered','trait.tsv'):return True
    if re.fullmatch(r'diagnostic_(?:8000|12000)_(?:rust|r)',sub[0]):
        return name not in {'cor_probe.tsv','tom_probe.tsv'} and p.suffix in {'.csv','.tsv','.txt'}
    if name in {'timings.tsv','module_sizes.tsv','module_colors.txt','filtered_genes.txt','module_trait_cor.tsv','module_trait_p.tsv'}:
        return 'negative_controls' not in sub and not sub[0].startswith('matrix_')
    return False

def main():
    entries=[];paths=[]
    for base in [ROOT/'validation/audit',ROOT/'validation/fixes']:
        for p in sorted(base.rglob('*')):
            if p.is_symlink():
                entries.append(dict(path=p.relative_to(ROOT).as_posix(),kind='local-symlink',published=False));continue
            if not p.is_file() or '__pycache__' in p.parts:continue
            choose=selected(p)
            entries.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),published=choose))
            if choose:paths.append(p.relative_to(ROOT).as_posix())
    manifest=dict(date='2026-10-28',performed_by='OpenAI Codex (ChatGPT)',policy='Selected review evidence; raw matrices, redundant exports and machine-specific symlinks remain local. This catalog does not imply full-square-matrix validation.',files=entries)
    (ROOT/'validation/PUBLISHED_EVIDENCE.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    print('\n'.join(paths))
if __name__=='__main__':main()
