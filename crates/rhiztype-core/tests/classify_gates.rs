use rhiztype_core::{RhizTypeEngine, ScreeningConfig, SequenceResult};

fn hxb2() -> String {
    let t = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/data/profiles_9719.tsv")).unwrap();
    for l in t.lines() {
        let mut it = l.split('\t');
        if it.next().map(|n| n.contains("HXB2")).unwrap_or(false) {
            return it.next().unwrap().chars().filter(|c| *c != '-' && *c != '?').collect();
        }
    }
    panic!("no HXB2 row");
}

fn run(a: usize, b: usize, mask: bool) -> SequenceResult {
    let seq = hxb2();
    let s = &seq[a - 1..b];
    let eng = RhizTypeEngine::new();
    let cfg = ScreeningConfig { mask_hypervariable_loops: mask, ..Default::default() };
    eng.classify_sequence("t", s, &cfg)
}

#[test]
fn fully_masked_slice_is_insufficient() {
    let r = run(6615, 6812, true);
    assert_eq!(r.locus_scope, "insufficient_coverage", "{:?}", r.structural_formula);
    assert!(r.segments.is_empty());
}

#[test]
fn span_start_excludes_masked_loop() {
    let r = run(6700, 7600, true);
    assert!(r.locus_scope == "insufficient_coverage" || r.hxb2_span_start >= 6813, "start {}", r.hxb2_span_start);
}

#[test]
fn gag_only_amplicon_is_not_pol() {
    let r = run(1800, 2200, true);
    assert_ne!(r.locus_scope, "pol", "{}-{} {}", r.hxb2_span_start, r.hxb2_span_end, r.alignment_mode);
}

#[test]
fn unmasked_slice_still_classified() {
    let r = run(6615, 6812, false);
    assert_ne!(r.locus_scope, "insufficient_coverage");
}
