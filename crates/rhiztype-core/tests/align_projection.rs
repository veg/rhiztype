use rhiztype_core::{RhizTypeEngine, ScreeningConfig, SequenceResult};

fn load_profile(name: &str) -> String {
    let tsv = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/profiles_9719.tsv"));
    for line in tsv.lines() {
        let mut p = line.trim().split('\t');
        if p.next() == Some(name) {
            return p.next().unwrap().to_string();
        }
    }
    panic!("profile {name} not found");
}

fn ungap(s: &str) -> String {
    s.chars().filter(|&c| c != '-' && c != '?').collect()
}

fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            o => o,
        })
        .collect()
}

fn slice(s: &str, a: usize, b: usize) -> String {
    s[a..b].to_string()
}

fn lcg_seq(n: usize, mut state: u64) -> String {
    (0..n)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[((state >> 33) & 3) as usize] as char
        })
        .collect()
}

fn run(eng: &RhizTypeEngine, seq: &str) -> SequenceResult {
    eng.classify_sequence("q", seq, &ScreeningConfig::default())
}

fn near(a: usize, b: usize, tol: usize) -> bool {
    a.abs_diff(b) <= tol
}

#[test]
fn a_exact_substring_span_end() {
    let eng = RhizTypeEngine::new();
    let hxb2 = load_profile("HXB2");
    for &(a, b) in &[(4107usize, 5107usize), (2000, 2300), (6300, 6450), (8000, 9100)] {
        let r = run(&eng, &slice(&hxb2, a, b));
        assert!(
            near(r.hxb2_span_end, b, 2),
            "[{a}..{b}] span_end {} != {b} ({})",
            r.hxb2_span_end,
            r.alignment_mode
        );
        assert!(near(r.hxb2_span_start, a + 1, 2), "[{a}..{b}] start {}", r.hxb2_span_start);
    }
}

#[test]
fn b_reverse_complement_matches_forward() {
    let eng = RhizTypeEngine::new();
    let hxb2 = load_profile("HXB2");
    let fwd = slice(&hxb2, 2000, 3500);
    let rf = run(&eng, &fwd);
    let rr = run(&eng, &revcomp(&fwd));
    assert!(rr.alignment_mode.contains("reverse-complemented"), "{}", rr.alignment_mode);
    assert!(!rf.alignment_mode.contains("reverse-complemented"));
    assert!(near(rf.hxb2_span_start, rr.hxb2_span_start, 2));
    assert!(near(rf.hxb2_span_end, rr.hxb2_span_end, 2));
}

#[test]
fn c_random_sequence_is_unclassifiable() {
    let eng = RhizTypeEngine::new();
    for seed in [1u64, 2, 3] {
        let r = run(&eng, &lcg_seq(1500, seed));
        assert_eq!(r.locus_scope, "insufficient_coverage", "seed {seed}: {}", r.structural_formula);
        assert!(r.matches.is_empty());
    }
}

#[test]
fn d_determinism() {
    let eng = RhizTypeEngine::new();
    let hxb2 = load_profile("HXB2");
    let rc = revcomp(&slice(&hxb2, 2000, 3500));
    let demos = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/examples/demos.fasta"));
    let recs = RhizTypeEngine::parse_fasta(demos);
    let cpx = recs.iter().find(|(h, _)| h.contains("06_cpx")).expect("CRF06_cpx").1.clone();
    for q in [rc, cpx] {
        let first = serde_json::to_string(&run(&eng, &q)).unwrap();
        for _ in 0..15 {
            assert_eq!(first, serde_json::to_string(&run(&eng, &q)).unwrap());
        }
    }
}

#[test]
fn e_rotated_row_is_not_prealigned() {
    let eng = RhizTypeEngine::new();
    let hxb2 = load_profile("HXB2");
    assert_eq!(hxb2.len(), 9719);
    let rot = format!("{}{}", &hxb2[30..], &hxb2[..30]);
    let r = run(&eng, &rot);
    assert!(!r.alignment_mode.starts_with("Pre-aligned"), "{}", r.alignment_mode);
    assert!(near(r.hxb2_span_start, 790, 3), "start {}", r.hxb2_span_start);
    let r0 = run(&eng, &hxb2);
    assert_eq!(r0.alignment_mode, "Pre-aligned (HXB2 9,719 nt)");
    assert_eq!(r0.structural_formula, "B(Clonal)");
}

#[test]
fn f_deletion_bases_land_on_true_columns() {
    let eng = RhizTypeEngine::new();
    let hxb2 = load_profile("HXB2");
    let q = format!("{}{}", &hxb2[900..1100], &hxb2[1130..1930]);
    assert_eq!(q.len(), 1000);
    let (cod, mode, _) = eng.project_to_hxb2_coding(&q);
    assert!(mode.starts_with("Auto-aligned"), "{mode}");
    let hb = hxb2.as_bytes();
    let (mut ok, mut tot) = (0, 0);
    for p0 in 1130..1930usize {
        let idx = p0 - 789;
        tot += 1;
        let b = match hb[p0] {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            _ => 3,
        };
        if cod[idx] == b {
            ok += 1;
        }
    }
    assert!(ok as f64 / tot as f64 > 0.99, "{ok}/{tot}");
}

#[test]
fn g_overlong_input_rejected() {
    let eng = RhizTypeEngine::new();
    let big = lcg_seq(20500, 7);
    let (cod, mode, _) = eng.project_to_hxb2_coding(&big);
    assert!(cod.iter().all(|&b| b < 0));
    assert!(mode.contains("rejected: >20000 nt"), "{mode}");
}

#[test]
fn h_tail_truncation_flag() {
    let eng = RhizTypeEngine::new();
    let hxb2 = load_profile("HXB2");
    let q = slice(&hxb2, 1000, 4000);
    let _ = ungap(&q);
    let r = run(&eng, &q);
    // exact substring: anchors reach the end, so no truncation flag expected
    assert!(!r.alignment_mode.contains("truncated"), "{}", r.alignment_mode);
}

#[test]
fn i_long_unanchored_tail_is_flagged() {
    let eng = RhizTypeEngine::new();
    let hxb2 = load_profile("HXB2");
    let q = format!("{}{}", &hxb2[1000..3000], lcg_seq(800, 11));
    let (_, mode, _) = eng.project_to_hxb2_coding(&q);
    assert!(mode.contains("3' tail truncated after 400 nt"), "{mode}");
}
