use rhiztype_core::*;

const PROFILES: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/profiles_9719.tsv"));

fn load_profile(name: &str) -> String {
    PROFILES
        .lines()
        .find_map(|l| {
            let (n, s) = l.split_once('\t')?;
            (n == name).then(|| s.trim().to_string())
        })
        .unwrap()
}

fn example(name: &str) -> String {
    std::fs::read_to_string(format!(
        "{}/../../docs/examples/{}.fasta",
        env!("CARGO_MANIFEST_DIR"),
        name
    ))
    .unwrap()
}

fn barcode_segs(bc: &str) -> Vec<(String, usize, usize)> {
    bc.split('|')
        .filter_map(|p| {
            let (t, c) = p.split_once(':')?;
            let c = c.split('[').next()?;
            let (s, e) = c.split_once('-')?;
            Some((t.to_string(), s.parse().ok()?, e.parse().ok()?))
        })
        .collect()
}

fn barcode_bps(bc: &str) -> Vec<(String, String, usize)> {
    let s = barcode_segs(bc);
    (0..s.len().saturating_sub(1))
        .map(|i| (s[i].0.clone(), s[i + 1].0.clone(), s[i + 1].1))
        .collect()
}

fn strict_shared(r: &SequenceResult, bc: &str) -> usize {
    let eb = barcode_bps(bc);
    r.breakpoints
        .iter()
        .filter(|b| {
            eb.iter().any(|(l, rt, p)| {
                *l == b.left_trunk && *rt == b.right_trunk && p.abs_diff(b.bp_hxb2) <= 200
            })
        })
        .count()
}

#[test]
fn gene_labels_gag_ca_and_pr_in() {
    assert_eq!(hxb2_to_short_gene(1171, 1893), "gag.CA");
    assert_eq!(hxb2_to_short_gene(2530, 4245), "PR..IN");
    assert_eq!(hxb2_to_short_gene(2290, 5096), "pol");
}

#[test]
fn compendium_donor_labels_survive_in_formulas() {
    let eng = RhizTypeEngine::new();
    for e in &eng.compendium {
        for (t, s, en) in barcode_segs(&e.hxb2_barcode) {
            if t == e.backbone || t == "?E" && false {
                continue;
            }
            let lab = hxb2_to_short_gene(s, en);
            if e.structural_formula.contains(&format!("{}:", t)) {
                assert!(
                    e.structural_formula.contains(&lab),
                    "{}: label {} for {}:{}-{} not in {}",
                    e.label, lab, t, s, en, e.structural_formula
                );
            }
        }
    }
}

#[test]
fn breakpoint_counts_are_strict_and_bounded() {
    let eng = RhizTypeEngine::new();
    let cfg = ScreeningConfig::default();
    let mut seen = 0;
    for name in ["demos", "crf01_02", "partial_env", "subtype_g"] {
        for (h, s) in RhizTypeEngine::parse_fasta(&example(name)) {
            let r = eng.classify_sequence(&h, &s, &cfg);
            for m in &r.matches {
                if let Some(rest) = m.match_type.split(", ").nth(1) {
                    if let Some(km) = rest.strip_suffix(" breakpoints)") {
                        let (k, mm) = km.split_once('/').unwrap();
                        let (k, mm): (usize, usize) = (k.parse().unwrap(), mm.parse().unwrap());
                        assert!(k <= mm, "{} {}: {}", h, m.label, m.match_type);
                        assert_eq!(k, strict_shared(&r, &m.hxb2_barcode), "{} {}: {}", h, m.label, m.match_type);
                        seen += 1;
                    }
                }
                if m.match_type.starts_with("Homologous Breakpoints") {
                    assert_eq!(strict_shared(&r, &m.hxb2_barcode), r.n_breakpoints, "{} {}", h, m.label);
                }
            }
        }
    }
    assert!(seen > 0);
}

#[test]
fn clonal_query_never_gets_mosaic_fallback() {
    let eng = RhizTypeEngine::new();
    let cm = load_profile("CM240");
    let cm: Vec<u8> = cm.bytes().collect();
    let mut q = vec![b'-'; 8628];
    for hx in 5097..=9417usize {
        q[hx - 790] = cm[hx - 1];
    }
    let q = String::from_utf8(q).unwrap();
    let r = eng.classify_sequence("cm240_partial", &q, &ScreeningConfig::default());
    if r.n_breakpoints == 0 {
        assert_ne!(r.matches[0].label, "Novel_URF_or_Unindexed_Mosaic", "{:?}", r.matches);
    }
}
