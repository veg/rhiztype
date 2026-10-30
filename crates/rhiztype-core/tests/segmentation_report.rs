use rhiztype_core::{RhizTypeEngine, ScreeningConfig, SequenceResult};

fn load_profile(name: &str) -> String {
    let txt = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/data/profiles_9719.tsv")).unwrap();
    for line in txt.lines() {
        let mut it = line.split('\t');
        if it.next() == Some(name) {
            return it.next().unwrap().to_string();
        }
    }
    panic!("profile {name} not found");
}

fn check_invariants(r: &SequenceResult) {
    for b in &r.breakpoints {
        assert_eq!(b.ci_width_nt, b.ci_hxb2_high - b.ci_hxb2_low + 1, "{} ci width", r.sequence_id);
    }
    for w in r.segments.windows(2) {
        assert_ne!(w[0].trunk, w[1].trunk, "{} adjacent same label", r.sequence_id);
    }
    assert_eq!(r.is_clonal, r.n_segments == 1, "{} is_clonal", r.sequence_id);
    if r.constituent_trunks.len() == 1 {
        assert_eq!(r.n_segments, 1, "{} single trunk but multiple segments", r.sequence_id);
    }
    if r.n_segments > 0 {
        assert_eq!(r.n_breakpoints, r.n_segments - 1);
    }
}

#[test]
fn examples_invariants() {
    let engine = RhizTypeEngine::new();
    let cfg = ScreeningConfig::default();
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/examples");
    let mut n = 0;
    for ent in std::fs::read_dir(dir).unwrap() {
        let p = ent.unwrap().path();
        if p.extension().and_then(|e| e.to_str()) != Some("fasta") {
            continue;
        }
        let txt = std::fs::read_to_string(&p).unwrap();
        for (h, s) in RhizTypeEngine::parse_fasta(&txt) {
            check_invariants(&engine.classify_sequence(&h, &s, &cfg));
            n += 1;
        }
    }
    assert!(n > 0);
}

#[test]
fn prealigned_hxb2_cm240_chimera_invariants() {
    let engine = RhizTypeEngine::new();
    let hxb2 = load_profile("HXB2");
    let cm = load_profile("CM240");
    let mut q = String::new();
    q.push_str(&hxb2[..5096]);
    q.push_str(&cm[5096..9417]);
    q.push_str(&hxb2[9417..]);
    assert_eq!(q.len(), 9719);
    let r = engine.classify_sequence("chimera", &q, &ScreeningConfig::default());
    check_invariants(&r);
}
