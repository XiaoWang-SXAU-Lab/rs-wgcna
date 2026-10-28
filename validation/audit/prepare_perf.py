from pathlib import Path
root=Path(__file__).resolve().parents[2]
s=(root/'examples/real_bench.rs').read_text();old='''let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);''';new='''let threads: usize = std::env::var("AUDIT_THREADS").unwrap_or("8".into()).parse().unwrap();
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();'''
assert old in s;s=s.replace(old,new)
(root/'examples/audit_perf.rs').write_text('// TEST-ONLY benchmark with explicit thread control; otherwise identical timing boundaries.\n'+s)
