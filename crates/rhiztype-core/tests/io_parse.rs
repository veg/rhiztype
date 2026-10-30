use rhiztype_core::{RhizTypeEngine, ScreeningConfig};

#[test]
fn bom_single_record() {
    let r = RhizTypeEngine::parse_fasta("\u{FEFF}>a\nACGT\n");
    assert_eq!(r, vec![("a".to_string(), "ACGT".to_string())]);
}

#[test]
fn cr_only_line_endings() {
    let r = RhizTypeEngine::parse_fasta(">a\rAC\rGT\r>b\rTT\r");
    assert_eq!(
        r,
        vec![
            ("a".to_string(), "ACGT".to_string()),
            ("b".to_string(), "TT".to_string())
        ]
    );
}

#[test]
fn headerless_paste() {
    let r = RhizTypeEngine::parse_fasta("acgt\nacgt\n");
    assert_eq!(r, vec![("unnamed_1".to_string(), "ACGTACGT".to_string())]);
}

#[test]
fn headerless_then_header() {
    let r = RhizTypeEngine::parse_fasta("AC\n>b\nGG\n");
    assert_eq!(r.len(), 2);
    assert_eq!(r[0].0, "unnamed_1");
    assert_eq!(r[1], ("b".to_string(), "GG".to_string()));
}

#[test]
fn validate_config() {
    assert!(ScreeningConfig::default().validate().is_ok());
    let bad = [
        ScreeningConfig { gamma: Some(f64::NAN), ..Default::default() },
        ScreeningConfig { gamma: Some(0.0), ..Default::default() },
        ScreeningConfig { gamma: Some(-1.0), ..Default::default() },
        ScreeningConfig { gamma: Some(f64::INFINITY), ..Default::default() },
        ScreeningConfig { p_fisher_thresh: 0.0, ..Default::default() },
        ScreeningConfig { p_fisher_thresh: f64::NAN, ..Default::default() },
        ScreeningConfig { p_fisher_thresh: 1.5, ..Default::default() },
        ScreeningConfig { min_span_nt: Some(0), ..Default::default() },
    ];
    for c in &bad {
        assert!(c.validate().is_err(), "{:?}", c);
    }
    assert!(ScreeningConfig { p_fisher_thresh: 1.0, gamma: Some(3.0), min_span_nt: Some(10), ..Default::default() }.validate().is_ok());
}
