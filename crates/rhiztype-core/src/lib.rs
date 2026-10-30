//! Core HXB2-anchored HIV-1 Group M subtyping and mosaic deconvolution engine.
//!
//! Supports both pre-aligned (HXB2 9,719-nt, coding 8,628-nt, or LANL 2023 11,531-nt MSA)
//! and unaligned FASTA sequences (full-length genomes, partial `pol`, partial `env`, or arbitrary amplicons)
//! via seeded 15-mer diagonal anchor hashing + banded Needleman-Wunsch coordinate projection onto HXB2 (`K03455`).

use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap, HashSet};

const PROFILES_TSV: &str = include_str!("../data/profiles_9719.tsv");
const COMPENDIUM_JSON: &str = include_str!("../data/compendium_bundle.json");

pub const HXB2_FULL_LEN: usize = 9719;
pub const CODING_START_HXB2: usize = 790;
pub const CODING_END_HXB2: usize = 9417;
pub const CODING_LEN: usize = CODING_END_HXB2 - CODING_START_HXB2 + 1; // 8628

/// Hypervariable indel-rich loops in `gag` (p6 PTAP) and `env` (`V1/V2`, `V3`, `V4`, `V5`)
const HYPERVARIABLE_MASK_RANGES: &[(usize, usize)] = &[
    (2130, 2185),
    (6615, 6812),
    (7110, 7217),
    (7377, 7478),
    (7596, 7637),
];

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CompendiumEntry {
    pub label: String,
    pub short_code: String,
    pub category: String,
    pub status: String,
    pub descriptive_name: String,
    pub structural_formula: String,
    pub backbone: String,
    pub t_mrca: Option<f64>,
    pub t_mrca_ci_low: Option<f64>,
    pub t_mrca_ci_high: Option<f64>,
    pub hxb2_barcode: String,
    pub scope: String,
}

#[derive(Debug, Clone, Deserialize)]
struct CompendiumBundle {
    lanl_msa_len: usize,
    lanl_msa_to_hxb2: Vec<usize>,
    compendium: Vec<CompendiumEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScreeningConfig {
    /// Template-switch penalty gamma (default 7.0)
    pub gamma: Option<f64>,
    /// Minimum segment span in HXB2 nucleotides (default 120 for full, 150 for partial)
    pub min_span_nt: Option<usize>,
    /// Minimum subtype-informative SNPs per segment (default 5)
    pub min_informative_snps: usize,
    /// Fisher exact test p-value threshold (default 1e-4)
    pub p_fisher_thresh: f64,
    /// Mask hypervariable gp120 V1-V5 loops and gag p6 PTAP duplication zone
    pub mask_hypervariable_loops: bool,
}

impl ScreeningConfig {
    /// Reject parameter values that would make the screen meaningless or non-deterministic.
    pub fn validate(&self) -> Result<(), String> {
        if let Some(g) = self.gamma {
            if !g.is_finite() || g <= 0.0 {
                return Err(format!("gamma must be a finite number > 0 (got {})", g));
            }
        }
        if self.min_span_nt == Some(0) {
            return Err("min-span must be >= 1".to_string());
        }
        if !self.p_fisher_thresh.is_finite() || self.p_fisher_thresh <= 0.0 || self.p_fisher_thresh > 1.0 {
            return Err(format!(
                "p-thresh must be in (0, 1] (got {})",
                self.p_fisher_thresh
            ));
        }
        Ok(())
    }
}

impl Default for ScreeningConfig {
    fn default() -> Self {
        Self {
            gamma: None,
            min_span_nt: None,
            min_informative_snps: 5,
            p_fisher_thresh: 1e-4,
            mask_hypervariable_loops: true,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SegmentCall {
    pub seg_idx: usize,
    pub trunk: String,
    pub hxb2_start: usize,
    pub hxb2_end: usize,
    pub span_nt: usize,
    pub ci_hxb2_start: usize,
    pub ci_hxb2_end: usize,
    pub short_gene: String,
    pub identity_pct: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BreakpointCall {
    pub bp_hxb2: usize,
    pub ci_hxb2_low: usize,
    pub ci_hxb2_high: usize,
    pub ci_width_nt: usize,
    pub left_trunk: String,
    pub right_trunk: String,
    pub p_fisher: f64,
    pub n11: usize,
    pub n12: usize,
    pub n21: usize,
    pub n22: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageMatch {
    pub label: String,
    pub descriptive_name: String,
    pub status: String,
    pub match_type: String,
    pub hxb2_barcode: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceResult {
    pub sequence_id: String,
    pub raw_length_nt: usize,
    pub mapped_coding_nt: usize,
    pub hxb2_span_start: usize,
    pub hxb2_span_end: usize,
    pub locus_scope: String,
    pub alignment_mode: String,
    pub structural_formula: String,
    pub descriptive_name: String,
    pub backbone_trunk: String,
    pub constituent_trunks: Vec<String>,
    pub is_clonal: bool,
    pub n_segments: usize,
    pub n_breakpoints: usize,
    pub hxb2_barcode: String,
    pub ci_bracketed_barcode: String,
    pub segments: Vec<SegmentCall>,
    pub breakpoints: Vec<BreakpointCall>,
    pub matches: Vec<LineageMatch>,
}

pub struct RhizTypeEngine {
    pub trunk_names: Vec<String>,
    /// 14 trunks x 8628 coding sites (0=A, 1=C, 2=G, 3=T, -1=indel, -2=forbidden domain)
    trunk_mat: Vec<Vec<i8>>,
    /// Full 9719-nt profiles for each trunk + HXB2 for k-mer anchoring and banded alignment
    full_profiles_9719: Vec<Vec<i8>>,
    cm240_coding: Vec<i8>,
    /// 15-mer hash index -> list of (profile_idx, hxb2_0based_pos)
    kmer_index: HashMap<u32, Vec<(u8, u16)>>,
    lanl_msa_len: usize,
    lanl_msa_to_hxb2: Vec<usize>,
    pub compendium: Vec<CompendiumEntry>,
    loop_mask_coding: Vec<bool>,
}

fn base_to_int(b: u8) -> i8 {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => -1,
    }
}

fn seq_to_int_vec(seq: &str) -> Vec<i8> {
    seq.bytes().map(base_to_int).collect()
}

fn encode_kmer_15(slice: &[i8]) -> Option<u32> {
    if slice.len() != 15 {
        return None;
    }
    let mut val: u32 = 0;
    for &b in slice {
        if b < 0 || b > 3 {
            return None;
        }
        val = (val << 2) | (b as u32);
    }
    Some(val)
}

impl Default for RhizTypeEngine {
    fn default() -> Self {
        Self::new()
    }
}

impl RhizTypeEngine {
    pub fn new() -> Self {
        let mut raw_profiles: HashMap<String, Vec<i8>> = HashMap::new();
        for line in PROFILES_TSV.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let mut parts = line.split('\t');
            if let (Some(name), Some(seq)) = (parts.next(), parts.next()) {
                raw_profiles.insert(name.to_string(), seq_to_int_vec(seq));
            }
        }

        let trunk_order = vec![
            "A1", "A2", "A6", "B", "C", "D", "F1", "F2", "G0", "H", "J", "K", "L", "?E",
        ];
        let mut trunk_names = Vec::new();
        let mut trunk_mat = Vec::new();
        let mut full_profiles_9719 = Vec::new();

        if let Some(hxb2_full) = raw_profiles.get("HXB2") {
            full_profiles_9719.push(hxb2_full.clone());
        }

        for &t_name in &trunk_order {
            let full_arr = raw_profiles
                .get(t_name)
                .expect("Missing trunk profile in embedded TSV");
            full_profiles_9719.push(full_arr.clone());
            let mut cod_arr = full_arr[(CODING_START_HXB2 - 1)..CODING_END_HXB2].to_vec();
            if t_name == "?E" {
                // ?E only exists from HXB2 5097..9417 (vif..nef)
                for idx in 0..cod_arr.len() {
                    let hxb2_pos = CODING_START_HXB2 + idx;
                    if hxb2_pos < 5097 {
                        cod_arr[idx] = -2;
                    }
                }
            }
            trunk_names.push(t_name.to_string());
            trunk_mat.push(cod_arr);
        }

        let cm240_full = raw_profiles
            .get("CM240")
            .expect("Missing CM240 profile in embedded TSV");
        let cm240_coding = cm240_full[(CODING_START_HXB2 - 1)..CODING_END_HXB2].to_vec();

        // Build 13-mer / 15-mer index across all 15 full profiles (step 2 nt for compact fast lookup)
        let mut kmer_index: HashMap<u32, Vec<(u8, u16)>> = HashMap::new();
        for (p_idx, prof) in full_profiles_9719.iter().enumerate() {
            for pos in (0..(prof.len().saturating_sub(15))).step_by(2) {
                if let Some(code) = encode_kmer_15(&prof[pos..(pos + 15)]) {
                    let entry = kmer_index.entry(code).or_default();
                    if entry.len() < 8 {
                        entry.push((p_idx as u8, pos as u16));
                    }
                }
            }
        }

        let bundle: CompendiumBundle =
            serde_json::from_str(COMPENDIUM_JSON).expect("Invalid embedded compendium_bundle.json");

        let mut loop_mask_coding = vec![false; CODING_LEN];
        for &(lo, hi) in HYPERVARIABLE_MASK_RANGES {
            for h in lo..=hi {
                if h >= CODING_START_HXB2 && h <= CODING_END_HXB2 {
                    loop_mask_coding[h - CODING_START_HXB2] = true;
                }
            }
        }

        Self {
            trunk_names,
            trunk_mat,
            full_profiles_9719,
            cm240_coding,
            kmer_index,
            lanl_msa_len: bundle.lanl_msa_len,
            lanl_msa_to_hxb2: bundle.lanl_msa_to_hxb2,
            compendium: bundle.compendium,
            loop_mask_coding,
        }
    }

    /// Parse a FASTA string into `(header, sequence)` tuples.
    pub fn parse_fasta(fasta_text: &str) -> Vec<(String, String)> {
        let mut records = Vec::new();
        let mut cur_hdr: Option<String> = None;
        let mut cur_seq = String::new();

        for line in fasta_text.trim_start_matches('\u{FEFF}').split(['\n', '\r']) {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            if let Some(rest) = line.strip_prefix('>') {
                if let Some(h) = cur_hdr.take() {
                    records.push((h, cur_seq.clone()));
                }
                cur_hdr = Some(rest.trim().to_string());
                cur_seq.clear();
            } else {
                if cur_hdr.is_none() {
                    cur_hdr = Some("unnamed_1".to_string());
                }
                for ch in line.chars() {
                    if !ch.is_whitespace() {
                        cur_seq.push(ch.to_ascii_uppercase());
                    }
                }
            }
        }
        if let Some(h) = cur_hdr {
            records.push((h, cur_seq));
        }
        records
    }

    /// Project any aligned or unaligned sequence onto the 8,628-nt HXB2 coding coordinate grid (790..9417).
    /// Returns `(coding_array, alignment_mode_string, raw_ungapped_len)`.
    pub fn project_to_hxb2_coding(&self, raw_seq: &str) -> (Vec<i8>, String, usize) {
        let seq_bytes = raw_seq.as_bytes();
        let total_len = seq_bytes.len();
        let ungapped: Vec<i8> = seq_bytes
            .iter()
            .filter(|&&b| b != b'-' && b != b'.')
            .map(|&b| base_to_int(b))
            .collect();
        let raw_ungapped_len = ungapped.len();
        let mut notes: Vec<String> = Vec::new();

        let layout_case = if total_len == HXB2_FULL_LEN {
            Some(1)
        } else if total_len == CODING_LEN {
            Some(2)
        } else if total_len == self.lanl_msa_len {
            Some(3)
        } else {
            None
        };

        if let Some(case) = layout_case {
            if self.prealigned_plausible(seq_bytes, case) {
                match case {
                    // Case 1: Exact 9,719-nt HXB2-anchored alignment
                    1 => {
                        let full_arr = seq_to_int_vec(raw_seq);
                        let cod = full_arr[(CODING_START_HXB2 - 1)..CODING_END_HXB2].to_vec();
                        return (cod, "Pre-aligned (HXB2 9,719 nt)".to_string(), raw_ungapped_len);
                    }
                    // Case 2: Exact 8,628-nt HXB2 coding alignment (790..9417)
                    2 => {
                        let cod = seq_to_int_vec(raw_seq);
                        return (
                            cod,
                            "Pre-aligned (HXB2 coding 8,628 nt)".to_string(),
                            raw_ungapped_len,
                        );
                    }
                    // Case 3: Exact 11,531-nt LANL 2023 Reference Alignment
                    _ => {
                        let full_arr = seq_to_int_vec(raw_seq);
                        let mut cod = vec![-1i8; CODING_LEN];
                        for (msa_col, &h_pos) in self.lanl_msa_to_hxb2.iter().enumerate() {
                            if h_pos >= CODING_START_HXB2 && h_pos <= CODING_END_HXB2 {
                                let c_idx = h_pos - CODING_START_HXB2;
                                if cod[c_idx] < 0 && full_arr[msa_col] >= 0 {
                                    cod[c_idx] = full_arr[msa_col];
                                }
                            }
                        }
                        return (
                            cod,
                            "Pre-aligned (LANL 2023 MSA 11,531 nt)".to_string(),
                            raw_ungapped_len,
                        );
                    }
                }
            }
            notes.push("length matched pre-aligned layout but anchors disagreed".to_string());
        }

        // Case 4: Unaligned or differently-aligned sequence -> Seed-and-extend coordinate projection onto HXB2
        let (full_9719, mut flags) = self.align_unaligned_to_hxb2_9719(&ungapped);
        notes.append(&mut flags);
        let cod = full_9719[(CODING_START_HXB2 - 1)..CODING_END_HXB2].to_vec();
        let mut mode = "Auto-aligned (15-mer anchor + banded NW to HXB2".to_string();
        for n in &notes {
            mode.push_str("; ");
            mode.push_str(n);
        }
        mode.push(')');
        (cod, mode, raw_ungapped_len)
    }

    /// Sanity check for records whose length equals a pre-aligned layout: gapped records are trusted;
    /// gap-free ones are probed with 15-mers against the HXB2/trunk index at the positions the layout predicts.
    fn prealigned_plausible(&self, seq: &[u8], case: u8) -> bool {
        if seq.iter().any(|&b| b == b'-' || b == b'.') {
            return true;
        }
        // (first column, last column exclusive) to probe, and the expected 0-based HXB2 position of a column
        let (lo, hi) = match case {
            1 => (CODING_START_HXB2 - 1, CODING_END_HXB2),
            2 => (0, CODING_LEN),
            _ => (0, seq.len()),
        };
        let expected = |c: usize| -> Option<usize> {
            match case {
                1 => Some(c),
                2 => Some(c + CODING_START_HXB2 - 1),
                _ => {
                    let h = *self.lanl_msa_to_hxb2.get(c)?;
                    if h >= CODING_START_HXB2 && h <= CODING_END_HXB2 {
                        Some(h - 1)
                    } else {
                        None
                    }
                }
            }
        };
        const N_PROBES: usize = 16;
        let span = hi.saturating_sub(lo).saturating_sub(40);
        let (mut agree, mut disagree) = (0usize, 0usize);
        for i in 0..N_PROBES {
            let start = lo + i * span / (N_PROBES - 1);
            // the index holds even reference positions only
            let mut c = start;
            let target = loop {
                if c >= start + 8 || c + 15 > seq.len() {
                    break None;
                }
                match expected(c) {
                    Some(e) if e % 2 == 0 => break Some(e),
                    _ => c += 1,
                }
            };
            let Some(target) = target else { continue };
            let bases: Vec<i8> = seq[c..c + 15].iter().map(|&b| base_to_int(b)).collect();
            let Some(code) = encode_kmer_15(&bases) else { continue };
            if let Some(locs) = self.kmer_index.get(&code) {
                if locs.iter().any(|&(_, r)| (r as usize).abs_diff(target) <= 3) {
                    agree += 1;
                } else if !locs.is_empty() {
                    disagree += 1;
                }
            }
        }
        disagree <= agree
    }

    /// Collect seed hits for one strand and chain them into non-overlapping anchors `(r_pos, q_pos, k)`.
    /// Returns the anchors and the chain score (`-1` when no chain exists).
    fn seed_and_chain(&self, q: &[i8]) -> (Vec<(usize, usize, usize)>, i32) {
        let qlen = q.len();
        // 1. Collect seed hits: (r_pos, q_pos, prof_idx, k)
        let mut hits: Vec<(usize, usize, usize, usize)> = Vec::new();
        for q_pos in 0..=(qlen.saturating_sub(15)) {
            if let Some(code) = encode_kmer_15(&q[q_pos..(q_pos + 15)]) {
                if let Some(locs) = self.kmer_index.get(&code) {
                    for &(p_idx, r_pos) in locs {
                        hits.push((r_pos as usize, q_pos, p_idx as usize, 15));
                    }
                }
            }
        }

        // If very few exact 15-mer hits, fall back to 11-mer scanning against the first 8 profiles
        // (HXB2, A1, A2, A6, B, C, D, F1)
        if hits.len() < 4 && qlen >= 80 {
            for (p_idx, prof) in self.full_profiles_9719.iter().enumerate().take(8) {
                let mut k11_map: HashMap<u32, usize> = HashMap::new();
                for r_pos in (0..(prof.len().saturating_sub(11))).step_by(3) {
                    let slice = &prof[r_pos..(r_pos + 11)];
                    if slice.iter().all(|&b| b >= 0 && b <= 3) {
                        let mut v = 0u32;
                        for &b in slice {
                            v = (v << 2) | (b as u32);
                        }
                        k11_map.entry(v).or_insert(r_pos);
                    }
                }
                for q_pos in (0..=(qlen.saturating_sub(11))).step_by(2) {
                    let slice = &q[q_pos..(q_pos + 11)];
                    if slice.iter().all(|&b| b >= 0 && b <= 3) {
                        let mut v = 0u32;
                        for &b in slice {
                            v = (v << 2) | (b as u32);
                        }
                        if let Some(&r_pos) = k11_map.get(&v) {
                            hits.push((r_pos, q_pos, p_idx, 11));
                        }
                    }
                }
            }
        }

        if hits.is_empty() {
            return (Vec::new(), -1);
        }

        // Deduplicate (r_pos, q_pos) keeping the lowest profile index (stable, so 15-mers win over 11-mers)
        hits.sort_by_key(|&(r, q_p, p_idx, _)| (q_p, r, p_idx));
        hits.dedup_by_key(|h| (h.0, h.1));

        // 2. Find dominant diagonal bin (r_pos - q_pos) to filter spurious off-diagonal repeats
        let mut diag_counts: BTreeMap<i32, usize> = BTreeMap::new();
        for &(r, qp, _, _) in &hits {
            let bin = ((r as i32) - (qp as i32)).div_euclid(250);
            *diag_counts.entry(bin).or_insert(0) += 1;
            *diag_counts.entry(bin - 1).or_insert(0) += 1;
            *diag_counts.entry(bin + 1).or_insert(0) += 1;
        }
        let best_bin = diag_counts
            .into_iter()
            .max_by_key(|&(bin, count)| (count, Reverse(bin.abs()), Reverse(bin)))
            .map(|(b, _)| b)
            .unwrap_or(0);
        let center_diag = best_bin * 250;

        let filtered: Vec<(usize, usize, usize)> = hits
            .into_iter()
            .filter(|&(r, qp, _, _)| ((r as i32 - qp as i32) - center_diag).abs() <= 600)
            .map(|(r, qp, _, k)| (r, qp, k))
            .collect();

        if filtered.is_empty() {
            return (Vec::new(), -1);
        }

        // 3. Collinear sparse DP chaining over anchors
        let n_hits = filtered.len();
        let mut dp: Vec<i32> = filtered.iter().map(|&(_, _, k)| k as i32).collect();
        let mut prev = vec![usize::MAX; n_hits];

        for i in 0..n_hits {
            let (ri, qi, ki) = filtered[i];
            let ki = ki as i32;
            let start_j = i.saturating_sub(120);
            for j in start_j..i {
                let (rj, qj, _) = filtered[j];
                if rj < ri && qj < qi {
                    let dr = (ri - rj) as i32;
                    let dq = (qi - qj) as i32;
                    let indel = (dr - dq).abs();
                    if indel <= 180 && dr <= 2500 && dq <= 2500 {
                        let overlap = (ki - dr.min(dq)).max(0);
                        let gain = ki - overlap - (indel / 3);
                        let cand = dp[j] + gain;
                        if cand > dp[i] {
                            dp[i] = cand;
                            prev[i] = j;
                        }
                    }
                }
            }
        }

        let mut best_end = 0;
        for i in 1..n_hits {
            if dp[i] > dp[best_end] {
                best_end = i;
            }
        }
        let score = dp[best_end];

        let mut chain: Vec<(usize, usize, usize)> = Vec::new();
        let mut cur = best_end;
        while cur != usize::MAX {
            chain.push(filtered[cur]);
            cur = prev[cur];
        }
        chain.reverse();

        // Merge overlapping anchors into non-overlapping anchor points
        let mut clean_anchors: Vec<(usize, usize, usize)> = Vec::new();
        for (r, qp, k) in chain {
            if let Some(&(lr, lq, lk)) = clean_anchors.last() {
                if r >= lr + lk && qp >= lq + lk {
                    clean_anchors.push((r, qp, k));
                }
            } else {
                clean_anchors.push((r, qp, k));
            }
        }
        (clean_anchors, score)
    }

    /// Aligns an unaligned sequence `q` onto the 9,719-nt HXB2 coordinate grid using 15-mer seed chaining
    /// and banded dynamic programming between collinear anchors. Both strands are chained and the
    /// stronger one is used. The second return value lists diagnostic notes for the mode string.
    fn align_unaligned_to_hxb2_9719(&self, q: &[i8]) -> (Vec<i8>, Vec<String>) {
        let qlen = q.len();
        let mut notes: Vec<String> = Vec::new();
        if qlen < 40 {
            return (vec![-1i8; HXB2_FULL_LEN], notes);
        }
        if qlen > 20000 {
            notes.push("rejected: >20000 nt".to_string());
            return (vec![-1i8; HXB2_FULL_LEN], notes);
        }

        let (fwd_anchors, fwd_score) = self.seed_and_chain(q);
        let rc: Vec<i8> = q
            .iter()
            .rev()
            .map(|&b| if (0..=3).contains(&b) { 3 - b } else { -1 })
            .collect();
        let (rc_anchors, rc_score) = self.seed_and_chain(&rc);

        let (clean_anchors, q, reverse) = if rc_score > fwd_score {
            (rc_anchors, &rc[..], true)
        } else {
            (fwd_anchors, q, false)
        };
        if clean_anchors.len() < 3 {
            return (vec![-1i8; HXB2_FULL_LEN], notes);
        }
        if reverse {
            notes.push("reverse-complemented".to_string());
        }

        let mut aligned = vec![-1i8; HXB2_FULL_LEN];
        let mut proportional_len = 0usize;

        // Fill left unanchored tail if short/moderate
        let (first_r, first_q, _) = clean_anchors[0];
        let mut left_untouched = 0usize;
        if first_q > 0 && first_r > 0 {
            let take_q = first_q.min(400);
            let take_r = first_r.min(take_q + 45);
            left_untouched = first_q - take_q;
            proportional_len = proportional_len.max(self.banded_align_segment(
                &q[(first_q - take_q)..first_q],
                first_r - take_r,
                first_r,
                true,
                false,
                &mut aligned,
            ));
        }

        // Fill exact anchors and banded gaps between consecutive anchors
        for idx in 0..clean_anchors.len() {
            let (r, qp, k) = clean_anchors[idx];
            for off in 0..k {
                if r + off < HXB2_FULL_LEN && qp + off < qlen {
                    aligned[r + off] = q[qp + off];
                }
            }
            if idx + 1 < clean_anchors.len() {
                let (next_r, next_q, _) = clean_anchors[idx + 1];
                let gap_r_start = r + k;
                let gap_q_start = qp + k;
                if next_r > gap_r_start && next_q > gap_q_start {
                    proportional_len = proportional_len.max(self.banded_align_segment(
                        &q[gap_q_start..next_q],
                        gap_r_start,
                        next_r,
                        false,
                        false,
                        &mut aligned,
                    ));
                }
            }
        }

        // Fill right unanchored tail
        let &(last_r, last_q, last_k) = clean_anchors.last().unwrap();
        let rem_q_start = last_q + last_k;
        let rem_r_start = last_r + last_k;
        let mut right_untouched = 0usize;
        if rem_q_start < qlen && rem_r_start < HXB2_FULL_LEN {
            let take_q = (qlen - rem_q_start).min(400);
            let take_r = (HXB2_FULL_LEN - rem_r_start).min(take_q + 45);
            right_untouched = qlen - rem_q_start - take_q;
            proportional_len = proportional_len.max(self.banded_align_segment(
                &q[rem_q_start..(rem_q_start + take_q)],
                rem_r_start,
                rem_r_start + take_r,
                false,
                true,
                &mut aligned,
            ));
        }

        if proportional_len > 0 {
            notes.push(format!("proportional projection over {} nt", proportional_len));
        }
        if left_untouched > 0 {
            notes.push(format!(
                "5' tail truncated after 400 nt ({} nt not aligned)",
                left_untouched
            ));
        }
        if right_untouched > 0 {
            notes.push(format!(
                "3' tail truncated after 400 nt ({} nt not aligned)",
                right_untouched
            ));
        }

        (aligned, notes)
    }

    /// Banded Needleman-Wunsch projection of `q_sub` onto HXB2 reference interval `[r_start..r_end]`,
    /// scoring matches against the consensus profile ensemble so any Group M subtype aligns cleanly.
    /// `free_ref_start` / `free_ref_end` make the reference slack at that end free (semi-global), leaving
    /// unused reference columns at -1. Returns the interval length when the interval was too large for
    /// DP and was projected proportionally, else 0.
    fn banded_align_segment(
        &self,
        q_sub: &[i8],
        r_start: usize,
        r_end: usize,
        free_ref_start: bool,
        free_ref_end: bool,
        aligned_out: &mut [i8],
    ) -> usize {
        let m = q_sub.len();
        let n = r_end.saturating_sub(r_start);
        if m == 0 || n == 0 {
            return 0;
        }
        // If interval is huge, project proportionally to avoid quadratic memory
        if m > 1200 || n > 1200 {
            for j in 0..n {
                let qi = (j * m) / n;
                aligned_out[r_start + j] = q_sub[qi];
            }
            return n.max(m);
        }

        let diff = (m as isize - n as isize).unsigned_abs();
        let band = diff + 35;
        let neg_inf = -1_000_000i32;

        let cols = n + 1;
        let mut score = vec![neg_inf; (m + 1) * cols];
        let mut trace = vec![0u8; (m + 1) * cols]; // 1=diag, 2=up (insertion in q), 3=left (deletion in q)

        score[0] = 0;
        for i in 1..=m.min(band) {
            score[i * cols] = -(i as i32) * 3;
            trace[i * cols] = 2;
        }
        for j in 1..=n.min(band) {
            score[j] = if free_ref_start { 0 } else { -(j as i32) * 3 };
            trace[j] = 3;
        }
        if free_ref_start {
            for j in 1..=n {
                score[j] = 0;
                trace[j] = 3;
            }
        }

        for i in 1..=m {
            let center_j = (i * n) / m;
            let j_lo = center_j.saturating_sub(band).max(1);
            let j_hi = (center_j + band).min(n);
            let qb = q_sub[i - 1];

            for j in j_lo..=j_hi {
                let r_pos = r_start + j - 1;
                // Match if query base appears in HXB2 or any clean trunk at r_pos
                let mut matched = false;
                if qb >= 0 {
                    for prof in &self.full_profiles_9719 {
                        if prof[r_pos] == qb {
                            matched = true;
                            break;
                        }
                    }
                }
                let sub_score = if matched { 2 } else { -1 };

                let diag = score[(i - 1) * cols + (j - 1)] + sub_score;
                let up = score[(i - 1) * cols + j] - 3;
                let left = score[i * cols + (j - 1)] - 3;

                let idx = i * cols + j;
                if diag >= up && diag >= left {
                    score[idx] = diag;
                    trace[idx] = 1;
                } else if up >= left {
                    score[idx] = up;
                    trace[idx] = 2;
                } else {
                    score[idx] = left;
                    trace[idx] = 3;
                }
            }
        }

        let mut i = m;
        let mut j = n;
        if free_ref_end {
            // Traceback from the best-scoring cell of the last row; later reference columns stay unaligned
            let j_lo = n.saturating_sub(band).max(1);
            let mut best = neg_inf;
            for jj in j_lo..=n {
                let s = score[m * cols + jj];
                if s >= best {
                    best = s;
                    j = jj;
                }
            }
        }
        while i > 0 && j > 0 {
            let t = trace[i * cols + j];
            if t == 1 {
                aligned_out[r_start + j - 1] = q_sub[i - 1];
                i -= 1;
                j -= 1;
            } else if t == 2 {
                i -= 1;
            } else if t == 3 {
                aligned_out[r_start + j - 1] = -1;
                j -= 1;
            } else {
                // Fallback if outside band
                aligned_out[r_start + j - 1] = q_sub[i - 1];
                i -= 1;
                j -= 1;
            }
        }
        0
    }

    /// Screen a single sequence `(header, raw_seq)` and return its full `SequenceResult`.
    pub fn classify_sequence(
        &self,
        header: &str,
        raw_seq: &str,
        cfg: &ScreeningConfig,
    ) -> SequenceResult {
        let (mut q_coding, align_mode, raw_len) = self.project_to_hxb2_coding(raw_seq);

        let mapped_coding_nt = q_coding.iter().filter(|&&b| b >= 0).count();

        // Apply artifact shield for hypervariable loops if enabled (before the span is determined)
        if cfg.mask_hypervariable_loops {
            for i in 0..CODING_LEN {
                if self.loop_mask_coding[i] {
                    q_coding[i] = -1;
                }
            }
        }

        // Determine active coding span [active_start..active_end] where usable nucleotides exist
        let valid_positions: Vec<usize> = q_coding
            .iter()
            .enumerate()
            .filter_map(|(i, &b)| if b >= 0 { Some(i) } else { None })
            .collect();

        if valid_positions.len() < 100 {
            return SequenceResult {
                sequence_id: header.to_string(),
                raw_length_nt: raw_len,
                mapped_coding_nt,
                hxb2_span_start: 0,
                hxb2_span_end: 0,
                locus_scope: "insufficient_coverage".to_string(),
                alignment_mode: align_mode,
                structural_formula: if cfg.mask_hypervariable_loops {
                    "Unclassifiable (<100 usable coding nt after loop mask)".to_string()
                } else {
                    "Unclassifiable (<100 mapped coding nt)".to_string()
                },
                descriptive_name: "Unclassifiable".to_string(),
                backbone_trunk: "None".to_string(),
                constituent_trunks: vec![],
                is_clonal: false,
                n_segments: 0,
                n_breakpoints: 0,
                hxb2_barcode: String::new(),
                ci_bracketed_barcode: String::new(),
                segments: vec![],
                breakpoints: vec![],
                matches: vec![],
            };
        }

        let first_idx = *valid_positions.first().unwrap();
        let last_idx = *valid_positions.last().unwrap();
        let h_start = CODING_START_HXB2 + first_idx;
        let h_end = CODING_START_HXB2 + last_idx;

        // Determine locus scope (full genome vs partial pol vs partial env vs custom partial)
        let span_width = h_end - h_start + 1;
        let locus_scope = if span_width >= 6800 && mapped_coding_nt >= 6000 {
            "full".to_string()
        } else if h_start >= 1800 && h_end > 2292 && h_end <= 5400 {
            "pol".to_string()
        } else if h_start >= 5900 && h_end <= 9100 {
            "env".to_string()
        } else {
            "partial".to_string()
        };

        let gamma = cfg.gamma.unwrap_or(7.0);
        let min_span_nt = cfg.min_span_nt.unwrap_or(if locus_scope == "full" {
            120
        } else {
            150
        });

        // Slice to [first_idx..=last_idx] so leading/trailing unsequenced regions in partial genes don't dilute edges
        let sub_len = last_idx - first_idx + 1;
        let q_sub = &q_coding[first_idx..=last_idx];
        let col_to_hxb2: Vec<usize> = (h_start..=h_end).collect();

        let (segs, bps) = self.run_paff_viterbi(
            q_sub,
            first_idx,
            &col_to_hxb2,
            sub_len,
            gamma,
            min_span_nt,
            cfg.min_informative_snps,
            cfg.p_fisher_thresh,
        );

        let (desc_name, backbone, struct_formula) =
            format_descriptive_name(&segs, &locus_scope);
        let hxb2_barcode = segs
            .iter()
            .map(|s| format!("{}:{}-{}", s.trunk, s.hxb2_start, s.hxb2_end))
            .collect::<Vec<_>>()
            .join("|");
        let ci_bracketed_barcode = segs
            .iter()
            .map(|s| {
                format!(
                    "{}:{}-{}[CI:{}-{}]",
                    s.trunk, s.hxb2_start, s.hxb2_end, s.ci_hxb2_start, s.ci_hxb2_end
                )
            })
            .collect::<Vec<_>>()
            .join("|");

        let mut constituent_set: HashSet<String> = HashSet::new();
        for s in &segs {
            constituent_set.insert(s.trunk.clone());
        }
        let mut constituent_trunks: Vec<String> = constituent_set.into_iter().collect();
        constituent_trunks.sort();

        let matches = self.match_compendium(&struct_formula, &segs, &bps, &locus_scope);

        SequenceResult {
            sequence_id: header.to_string(),
            raw_length_nt: raw_len,
            mapped_coding_nt,
            hxb2_span_start: h_start,
            hxb2_span_end: h_end,
            locus_scope,
            alignment_mode: align_mode,
            structural_formula: struct_formula,
            descriptive_name: desc_name,
            backbone_trunk: backbone,
            constituent_trunks,
            is_clonal: segs.len() == 1,
            n_segments: segs.len(),
            n_breakpoints: bps.len(),
            hxb2_barcode,
            ci_bracketed_barcode,
            segments: segs,
            breakpoints: bps,
            matches,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn run_paff_viterbi(
        &self,
        q_sub: &[i8],
        offset_in_coding: usize,
        col_to_hxb2: &[usize],
        l_sub: usize,
        gamma: f64,
        min_hxb2_span_nt: usize,
        min_informative_snps: usize,
        p_fisher_thresh: f64,
    ) -> (Vec<SegmentCall>, Vec<BreakpointCall>) {
        let k_trunks = self.trunk_names.len(); // 14
        let k_total = k_trunks + 1;
        let ghost_idx = k_trunks;
        let mismatch_cost = 1.0f64;
        let ghost_cost = 0.58f64;

        let mut dp = vec![0.0f64; k_total * l_sub];
        let mut back = vec![0u16; k_total * l_sub];

        // Site 0 initialization
        for k in 0..k_trunks {
            dp[k * l_sub] = self.site_emission_cost(
                q_sub[0],
                self.trunk_mat[k][offset_in_coding],
                mismatch_cost,
            );
        }
        dp[ghost_idx * l_sub] = if q_sub[0] >= 0 { ghost_cost } else { 0.0 };

        for u in 1..l_sub {
            let mut best_prev_idx = 0usize;
            let mut best_prev_val = dp[u - 1];
            for k in 1..k_total {
                let v = dp[k * l_sub + (u - 1)];
                if v < best_prev_val {
                    best_prev_val = v;
                    best_prev_idx = k;
                }
            }
            let switch_cost = best_prev_val + gamma;
            let cod_pos = offset_in_coding + u;

            for k in 0..k_total {
                let stay_cost = dp[k * l_sub + (u - 1)];
                let emit = if k < k_trunks {
                    let base_c =
                        self.site_emission_cost(q_sub[u], self.trunk_mat[k][cod_pos], mismatch_cost);
                    // Require >3.5% SNP advantage before preferring the partial ?E (CRF01_AE) cassette over sampled basal trunks
                    if k + 1 == k_trunks && q_sub[u] >= 0 {
                        base_c + 0.035
                    } else {
                        base_c
                    }
                } else if q_sub[u] >= 0 {
                    ghost_cost
                } else {
                    0.0
                };
                if switch_cost < stay_cost {
                    dp[k * l_sub + u] = switch_cost + emit;
                    back[k * l_sub + u] = best_prev_idx as u16;
                } else {
                    dp[k * l_sub + u] = stay_cost + emit;
                    back[k * l_sub + u] = k as u16;
                }
            }
        }

        let mut path = vec![0usize; l_sub];
        let mut best_last = 0usize;
        let mut best_last_val = dp[l_sub - 1];
        for k in 1..k_total {
            let v = dp[k * l_sub + (l_sub - 1)];
            if v < best_last_val {
                best_last_val = v;
                best_last = k;
            }
        }
        path[l_sub - 1] = best_last;
        for u in (1..l_sub).rev() {
            path[u - 1] = back[path[u] * l_sub + u] as usize;
        }

        let mut segs: Vec<(usize, usize, usize)> = Vec::new();
        let mut cur_t = path[0];
        let mut start_u = 0usize;
        for (u, &t) in path.iter().enumerate().skip(1) {
            if t != cur_t {
                segs.push((start_u, u, cur_t));
                cur_t = t;
                start_u = u;
            }
        }
        segs.push((start_u, l_sub, cur_t));

        // Iterative statistical boundary pruning
        let mut changed = true;
        while changed && segs.len() > 1 {
            changed = false;
            let mut best_merge_idx: Option<usize> = None;
            let mut worst_score = -1.0f64;

            for i in 0..(segs.len() - 1) {
                let (s1, e1, t1) = segs[i];
                let (s2, e2, t2) = segs[i + 1];
                if t1 == t2 {
                    best_merge_idx = Some(i);
                    break;
                }
                let span1 = (col_to_hxb2[e1 - 1] + 1).saturating_sub(col_to_hxb2[s1]);
                let span2 = (col_to_hxb2[e2 - 1] + 1).saturating_sub(col_to_hxb2[s2]);
                let (pval, n11, _n12, _n21, n22, _, _) = self.test_boundary_significance(
                    q_sub,
                    offset_in_coding,
                    ghost_idx,
                    s1,
                    e1,
                    e2,
                    t1,
                    t2,
                    min_informative_snps,
                );
                if span1 < min_hxb2_span_nt
                    || span2 < min_hxb2_span_nt
                    || pval > p_fisher_thresh
                    || n11 < min_informative_snps
                    || n22 < min_informative_snps
                {
                    let score = pval
                        + if span1 < min_hxb2_span_nt || span2 < min_hxb2_span_nt {
                            1.5
                        } else {
                            0.0
                        };
                    if score > worst_score {
                        worst_score = score;
                        best_merge_idx = Some(i);
                    }
                }
            }

            if let Some(i) = best_merge_idx {
                let (s1, e1, t1) = segs[i];
                let (s2, e2, t2) = segs[i + 1];
                let span1 = (col_to_hxb2[e1 - 1] + 1).saturating_sub(col_to_hxb2[s1]);
                let span2 = (col_to_hxb2[e2 - 1] + 1).saturating_sub(col_to_hxb2[s2]);

                let valid1 = self.trunk_valid_fraction(t1, offset_in_coding + s1, offset_in_coding + e2);
                let valid2 = self.trunk_valid_fraction(t2, offset_in_coding + s1, offset_in_coding + e2);

                let keep_t = if valid1 < 0.7 && valid2 >= 0.7 {
                    t2
                } else if valid2 < 0.7 && valid1 >= 0.7 {
                    t1
                } else if span1 >= span2 {
                    t1
                } else {
                    t2
                };
                segs[i] = (s1, e2, keep_t);
                segs.remove(i + 1);
                changed = true;
            }
        }

        // Merge identical adjacent states
        let mut state_merged: Vec<(usize, usize, usize)> = Vec::new();
        for seg in segs {
            if let Some(last) = state_merged.last_mut() {
                if last.2 == seg.2 {
                    last.1 = seg.1;
                    continue;
                }
            }
            state_merged.push(seg);
        }

        // Label each piece, then merge adjacent pieces whose final label is equal
        // (e.g. trunk ?E next to a ghost segment relabelled ?E).
        let label_of = |s: usize, e: usize, t: usize| -> String {
            if t < k_trunks {
                return self.trunk_names[t].clone();
            }
            let h0 = col_to_hxb2[s];
            let mut valid_cnt = 0usize;
            let mut match_cnt = 0usize;
            for u in s..e {
                let qb = q_sub[u];
                let cb = self.cm240_coding[offset_in_coding + u];
                if qb >= 0 && cb >= 0 {
                    valid_cnt += 1;
                    if qb == cb {
                        match_cnt += 1;
                    }
                }
            }
            if valid_cnt > 50 && h0 >= 5000 && (match_cnt as f64 / valid_cnt as f64) >= 0.85 {
                "?E".to_string()
            } else {
                "Unattributed".to_string()
            }
        };
        // (start, end, state, label, mixed states, span of the piece keeping the state)
        let mut merged: Vec<(usize, usize, usize, String, bool, usize)> = Vec::new();
        for &(s, e, t) in &state_merged {
            let label = label_of(s, e, t);
            let span = (col_to_hxb2[e - 1] + 1).saturating_sub(col_to_hxb2[s]);
            if let Some(last) = merged.last_mut() {
                if last.3 == label {
                    last.4 |= last.2 != t;
                    if span > last.5 {
                        last.2 = t;
                        last.5 = span;
                    }
                    last.1 = e;
                    continue;
                }
            }
            merged.push((s, e, t, label, false, span));
        }

        let mut breakpoints = Vec::new();
        for i in 0..(merged.len().saturating_sub(1)) {
            let (s1, e1, t1, _, _, _) = merged[i];
            let (_s2, e2, t2, _, _, _) = merged[i + 1];
            let (pval, n11, n12, n21, n22, ci_low_col, ci_high_col) = self
                .test_boundary_significance(
                    q_sub,
                    offset_in_coding,
                    ghost_idx,
                    s1,
                    e1,
                    e2,
                    t1,
                    t2,
                    min_informative_snps,
                );
            let ci_low_h = col_to_hxb2[ci_low_col];
            let ci_high_h = col_to_hxb2[ci_high_col];
            breakpoints.push(BreakpointCall {
                bp_hxb2: col_to_hxb2[e1],
                ci_hxb2_low: ci_low_h,
                ci_hxb2_high: ci_high_h,
                ci_width_nt: ci_high_h.saturating_sub(ci_low_h) + 1,
                left_trunk: String::new(),
                right_trunk: String::new(),
                p_fisher: pval,
                n11,
                n12,
                n21,
                n22,
            });
        }

        let mut out_segs = Vec::new();
        for (i, (s, e, t, t_label, mixed, _)) in merged.iter().cloned().enumerate() {
            let h0 = col_to_hxb2[s];
            let h1 = col_to_hxb2[e - 1];
            let ci_start = if i > 0 {
                breakpoints[i - 1].ci_hxb2_low
            } else {
                h0
            };
            let ci_end = if i < breakpoints.len() {
                breakpoints[i].ci_hxb2_high
            } else {
                h1
            };

            // Compute identity % against assigned trunk
            let mut id_valid = 0usize;
            let mut id_match = 0usize;
            let t_check = if t < k_trunks && !mixed {
                Some(&self.trunk_mat[t])
            } else if t_label == "?E" {
                Some(&self.cm240_coding)
            } else {
                None
            };
            if let Some(ref_arr) = t_check {
                for u in s..e {
                    let qb = q_sub[u];
                    let rb = ref_arr[offset_in_coding + u];
                    if qb >= 0 && rb >= 0 {
                        id_valid += 1;
                        if qb == rb {
                            id_match += 1;
                        }
                    }
                }
            }
            let identity_pct = if id_valid > 0 {
                (100.0 * id_match as f64 / id_valid as f64 * 10.0).round() / 10.0
            } else {
                0.0
            };

            out_segs.push(SegmentCall {
                seg_idx: i + 1,
                trunk: t_label,
                hxb2_start: h0,
                hxb2_end: h1,
                span_nt: h1 - h0 + 1,
                ci_hxb2_start: ci_start,
                ci_hxb2_end: ci_end,
                short_gene: hxb2_to_short_gene(h0, h1),
                identity_pct,
            });
        }

        for (i, bp) in breakpoints.iter_mut().enumerate() {
            bp.left_trunk = out_segs[i].trunk.clone();
            bp.right_trunk = out_segs[i + 1].trunk.clone();
        }

        (out_segs, breakpoints)
    }

    #[inline]
    fn site_emission_cost(&self, q_base: i8, ref_base: i8, mismatch_cost: f64) -> f64 {
        if q_base < 0 {
            return 0.0;
        }
        if ref_base == -2 {
            return 2.0;
        }
        if ref_base == -1 {
            return 0.35;
        }
        if q_base != ref_base {
            mismatch_cost
        } else {
            0.0
        }
    }

    fn trunk_valid_fraction(&self, t_idx: usize, cod_start: usize, cod_end: usize) -> f64 {
        if t_idx >= self.trunk_mat.len() || cod_end <= cod_start {
            return 1.0;
        }
        let slice = &self.trunk_mat[t_idx][cod_start..cod_end];
        let valid = slice.iter().filter(|&&b| b >= 0).count();
        valid as f64 / slice.len() as f64
    }

    #[allow(clippy::too_many_arguments)]
    fn test_boundary_significance(
        &self,
        q_sub: &[i8],
        offset: usize,
        ghost_idx: usize,
        col_start: usize,
        col_bp: usize,
        col_end: usize,
        t_l: usize,
        t_r: usize,
        min_snps: usize,
    ) -> (f64, usize, usize, usize, usize, usize, usize) {
        let mut n11 = 0usize;
        let mut n12 = 0usize;
        let mut n21 = 0usize;
        let mut n22 = 0usize;

        if t_l != ghost_idx && t_r != ghost_idx {
            let ref_l = &self.trunk_mat[t_l];
            let ref_r = &self.trunk_mat[t_r];

            let valid_r_on_left = (col_start..col_bp)
                .filter(|&u| ref_r[offset + u] >= 0)
                .count() as f64
                / ((col_bp - col_start).max(1) as f64);

            if valid_r_on_left < 0.5 {
                for u in col_start..col_bp {
                    let qb = q_sub[u];
                    let rl = ref_l[offset + u];
                    if qb >= 0 && rl >= 0 {
                        let other_diff = self
                            .trunk_mat
                            .iter()
                            .any(|tr| tr[offset + u] >= 0 && tr[offset + u] != rl);
                        if other_diff {
                            if qb == rl {
                                n11 += 1;
                            } else {
                                n12 += 1;
                            }
                        }
                    }
                }
            } else {
                for u in col_start..col_bp {
                    let qb = q_sub[u];
                    let rl = ref_l[offset + u];
                    let rr = ref_r[offset + u];
                    if qb >= 0 && rl >= 0 && rr >= 0 && rl != rr {
                        if qb == rl {
                            n11 += 1;
                        } else if qb == rr {
                            n12 += 1;
                        }
                    }
                }
            }

            let valid_l_on_right = (col_bp..col_end)
                .filter(|&u| ref_l[offset + u] >= 0)
                .count() as f64
                / ((col_end - col_bp).max(1) as f64);

            if valid_l_on_right < 0.5 {
                for u in col_bp..col_end {
                    let qb = q_sub[u];
                    let rr = ref_r[offset + u];
                    if qb >= 0 && rr >= 0 {
                        let other_diff = self
                            .trunk_mat
                            .iter()
                            .any(|tr| tr[offset + u] >= 0 && tr[offset + u] != rr);
                        if other_diff {
                            if qb == rr {
                                n22 += 1;
                            } else {
                                n21 += 1;
                            }
                        }
                    }
                }
            } else {
                for u in col_bp..col_end {
                    let qb = q_sub[u];
                    let rl = ref_l[offset + u];
                    let rr = ref_r[offset + u];
                    if qb >= 0 && rl >= 0 && rr >= 0 && rl != rr {
                        if qb == rl {
                            n21 += 1;
                        } else if qb == rr {
                            n22 += 1;
                        }
                    }
                }
            }
        } else {
            let known_t = if t_l == ghost_idx { t_r } else { t_l };
            let ref_k = &self.trunk_mat[known_t];

            let mut match_k_l = 0usize;
            let mut orphan_l = 0usize;
            for u in col_start..col_bp {
                let qb = q_sub[u];
                let rk = ref_k[offset + u];
                if qb >= 0 && rk >= 0 {
                    if qb == rk {
                        match_k_l += 1;
                    } else if !self.trunk_mat.iter().any(|tr| tr[offset + u] == qb) {
                        orphan_l += 1;
                    }
                }
            }

            let mut match_k_r = 0usize;
            let mut orphan_r = 0usize;
            for u in col_bp..col_end {
                let qb = q_sub[u];
                let rk = ref_k[offset + u];
                if qb >= 0 && rk >= 0 {
                    if qb == rk {
                        match_k_r += 1;
                    } else if !self.trunk_mat.iter().any(|tr| tr[offset + u] == qb) {
                        orphan_r += 1;
                    }
                }
            }

            if t_l == ghost_idx {
                n11 = orphan_l;
                n12 = match_k_l;
                n21 = orphan_r;
                n22 = match_k_r;
            } else {
                n11 = match_k_l;
                n12 = orphan_l;
                n21 = match_k_r;
                n22 = orphan_r;
            }
        }

        let pval = if n11 >= min_snps && n22 >= min_snps && (n11 * n22 > n12 * n21) {
            fisher_exact_greater(n11, n12, n21, n22)
        } else {
            1.0
        };

        // Informative site CI brackets
        let mut left_info_last = col_bp;
        let mut right_info_first = col_bp;

        if t_l != ghost_idx && t_r != ghost_idx {
            let ref_l = &self.trunk_mat[t_l];
            let ref_r = &self.trunk_mat[t_r];
            for u in col_start..col_bp {
                let qb = q_sub[u];
                let rl = ref_l[offset + u];
                let rr = ref_r[offset + u];
                if qb >= 0 && rl >= 0 && rr >= 0 && rl != rr && (qb == rl || qb == rr) {
                    left_info_last = u;
                }
            }
            for u in col_bp..col_end {
                let qb = q_sub[u];
                let rl = ref_l[offset + u];
                let rr = ref_r[offset + u];
                if qb >= 0 && rl >= 0 && rr >= 0 && rl != rr && (qb == rl || qb == rr) {
                    right_info_first = u;
                    break;
                }
            }
        }

        (
            pval,
            n11,
            n12,
            n21,
            n22,
            left_info_last,
            right_info_first,
        )
    }

    /// Match a sequence's structural formula, segment coordinates, and breakpoints against the 200-lineage compendium.
    fn match_compendium(
        &self,
        struct_formula: &str,
        segs: &[SegmentCall],
        bps: &[BreakpointCall],
        locus_scope: &str,
    ) -> Vec<LineageMatch> {
        let mut matches = Vec::new();
        let base_formula = struct_formula
            .trim_end_matches("@pol")
            .trim_end_matches("@env")
            .trim_end_matches("@partial");

        // 1. Exact structural formula matches
        for entry in &self.compendium {
            if entry.structural_formula == base_formula
                || (base_formula == "G0(Clonal)" && entry.label == "Subtype_G")
            {
                if (locus_scope == "full" && entry.scope == "full")
                    || (locus_scope == entry.scope)
                    || segs.len() == 1
                {
                    matches.push(LineageMatch {
                        label: entry.label.clone(),
                        descriptive_name: entry.descriptive_name.clone(),
                        status: entry.status.clone(),
                        match_type: if entry.structural_formula == base_formula {
                            "Exact Structural Formula".to_string()
                        } else {
                            "Basal Trunk Match (G0)".to_string()
                        },
                        hxb2_barcode: entry.hxb2_barcode.clone(),
                    });
                }
            }
        }

        // 2. Check coordinate-level mosaic concordance and breakpoint proximity (+/- 200 nt) against compendium barcodes
        if !bps.is_empty() && matches.len() < 4 {
            let mut scored_candidates: Vec<(f64, usize, usize, usize, &CompendiumEntry)> = Vec::new();

            for entry in &self.compendium {
                if matches.iter().any(|m| m.label == entry.label) {
                    continue;
                }
                let entry_segs = parse_barcode_segments(&entry.hxb2_barcode);
                let entry_bps = parse_barcode_breakpoints(&entry.hxb2_barcode);
                if entry_segs.is_empty() || entry_bps.is_empty() {
                    continue;
                }

                // Compute coordinate overlap concordance across query segments
                let mut total_ovl = 0usize;
                let mut same_trunk_ovl = 0usize;
                for q_seg in segs {
                    for (e_trk, e_s, e_e) in &entry_segs {
                        let ovl_s = q_seg.hxb2_start.max(*e_s);
                        let ovl_e = q_seg.hxb2_end.min(*e_e);
                        if ovl_e >= ovl_s {
                            let len = ovl_e - ovl_s + 1;
                            total_ovl += len;
                            if &q_seg.trunk == e_trk {
                                same_trunk_ovl += len;
                            }
                        }
                    }
                }
                if total_ovl < 400 {
                    continue;
                }
                let concordance = same_trunk_ovl as f64 / total_ovl as f64;

                let q_min_h = segs.first().map(|s| s.hxb2_start).unwrap_or(790);
                let q_max_h = segs.last().map(|s| s.hxb2_end).unwrap_or(9417);
                let span_ebps = entry_bps
                    .iter()
                    .filter(|(_, _, pos)| *pos >= q_min_h.saturating_sub(150) && *pos <= q_max_h + 150)
                    .count()
                    .max(1);

                let mut shared = 0usize;
                let mut pos_shared = 0usize;
                for bp in bps {
                    if entry_bps.iter().any(|(l_tr, r_tr, pos)| {
                        l_tr == &bp.left_trunk
                            && r_tr == &bp.right_trunk
                            && (*pos as isize - bp.bp_hxb2 as isize).unsigned_abs() <= 200
                    }) {
                        shared += 1;
                    }
                    if entry_bps.iter().any(|(l_tr, r_tr, pos)| {
                        (l_tr == &bp.left_trunk || r_tr == &bp.right_trunk)
                            && (*pos as isize - bp.bp_hxb2 as isize).unsigned_abs() <= 200
                    }) {
                        pos_shared += 1;
                    }
                }

                if (shared == bps.len() && shared == span_ebps)
                    || (concordance >= 0.88 && shared >= 1)
                    || (concordance >= 0.91 && pos_shared >= 2)
                    || (shared >= 3 && shared * 4 >= span_ebps * 3)
                {
                    scored_candidates.push((concordance, shared, pos_shared, span_ebps, entry));
                }
            }

            scored_candidates.sort_by(|a, b| {
                b.0.partial_cmp(&a.0)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| b.1.cmp(&a.1))
                    .then_with(|| b.2.cmp(&a.2))
                    .then_with(|| a.4.label.cmp(&b.4.label))
            });

            for (concordance, shared, _pos_shared, span_ebps, entry) in scored_candidates.into_iter().take(4) {
                let match_type = if shared == bps.len() && shared == span_ebps {
                    format!(
                        "Homologous Breakpoints ({:.1}% coordinate concordance)",
                        concordance * 100.0
                    )
                } else {
                    format!(
                        "Homologous Mosaic ({:.1}% concordance, {}/{} breakpoints)",
                        concordance * 100.0,
                        shared,
                        span_ebps.max(shared)
                    )
                };
                matches.push(LineageMatch {
                    label: entry.label.clone(),
                    descriptive_name: entry.descriptive_name.clone(),
                    status: entry.status.clone(),
                    match_type,
                    hxb2_barcode: entry.hxb2_barcode.clone(),
                });
            }
        }

        if matches.is_empty() && bps.is_empty() {
            matches.push(LineageMatch {
                label: "Unindexed_Clonal_Lineage".to_string(),
                descriptive_name: struct_formula.to_string(),
                status: "Clonal call with no compendium entry for this trunk/scope".to_string(),
                match_type: "No Compendium Match".to_string(),
                hxb2_barcode: String::new(),
            });
        } else if matches.is_empty() {
            matches.push(LineageMatch {
                label: "Novel_URF_or_Unindexed_Mosaic".to_string(),
                descriptive_name: struct_formula.to_string(),
                status: "Unique or Unindexed Recombinant Architecture".to_string(),
                match_type: "De Novo".to_string(),
                hxb2_barcode: String::new(),
            });
        }

        matches.truncate(6);
        matches
    }
}

fn parse_barcode_segments(barcode: &str) -> Vec<(String, usize, usize)> {
    let mut parsed_segs: Vec<(String, usize, usize)> = Vec::new();
    for part in barcode.split('|') {
        if let Some((trk, coords)) = part.split_once(':') {
            let coords_clean = coords.split('[').next().unwrap_or(coords);
            if let Some((s_str, e_str)) = coords_clean.split_once('-') {
                if let (Ok(s), Ok(e)) = (s_str.parse::<usize>(), e_str.parse::<usize>()) {
                    parsed_segs.push((trk.to_string(), s, e));
                }
            }
        }
    }
    parsed_segs
}

fn parse_barcode_breakpoints(barcode: &str) -> Vec<(String, String, usize)> {
    let parsed_segs = parse_barcode_segments(barcode);
    let mut bps = Vec::new();
    for i in 0..(parsed_segs.len().saturating_sub(1)) {
        bps.push((
            parsed_segs[i].0.clone(),
            parsed_segs[i + 1].0.clone(),
            parsed_segs[i + 1].1,
        ));
    }
    bps
}

/// Map an HXB2 nucleotide interval `[h0, h1]` to a concise HIV-1 gene/domain span label.
pub fn hxb2_to_short_gene(h0: usize, h1: usize) -> String {
    let span = h1.saturating_sub(h0);
    let pt_gene = |pos: usize| -> &'static str {
        if pos < 790 {
            "5'LTR"
        } else if pos <= 2292 {
            if pos <= 1185 {
                "gag.MA"
            } else if pos <= 1878 {
                "gag.CA"
            } else {
                "gag.NC"
            }
        } else if pos <= 5096 {
            if pos <= 2549 {
                "PR"
            } else if pos <= 4229 {
                "RT"
            } else {
                "IN"
            }
        } else if pos <= 5559 {
            "vif"
        } else if pos <= 5850 {
            "vpr"
        } else if pos <= 6062 {
            "tat1"
        } else if pos <= 6225 {
            "vpu"
        } else if pos <= 8795 {
            if pos <= 7757 {
                "gp120"
            } else {
                "gp41"
            }
        } else if pos <= 9417 {
            "nef"
        } else {
            "3'LTR"
        }
    };

    let g0 = pt_gene((h0 + 15).min(h1));
    let g1 = pt_gene(h1.saturating_sub(15).max(h0));

    if g0 == g1 {
        if (g0 == "PR" || g0 == "RT" || g0 == "IN") && span > 1800 {
            return "pol".to_string();
        }
        if (g0 == "gp120" || g0 == "gp41") && span > 1800 {
            return "env".to_string();
        }
        return g0.to_string();
    }
    if g0.starts_with("gag.") && g1.starts_with("gag.") {
        if span > 900 {
            return "gag".to_string();
        }
        return format!("{}..{}", g0, g1);
    }
    if g0 == "PR" && g1 == "IN" {
        if span > 1800 {
            return "pol".to_string();
        }
        return "PR..IN".to_string();
    }
    if g0 == "gp120" && g1 == "gp41" {
        if span > 1800 {
            return "env".to_string();
        }
        return "gp120..gp41".to_string();
    }
    let g0_norm = if g0.starts_with("gag.") && !g1.starts_with("gag.") {
        "gag"
    } else {
        g0
    };
    let g1_norm = if g1.starts_with("gag.") && !g0.starts_with("gag.") {
        "gag"
    } else {
        g1
    };
    format!("{}..{}", g0_norm, g1_norm)
}

fn format_descriptive_name(
    segs: &[SegmentCall],
    locus_scope: &str,
) -> (String, String, String) {
    let mut trunk_spans: HashMap<String, usize> = HashMap::new();
    for s in segs {
        *trunk_spans.entry(s.trunk.clone()).or_insert(0) += s.span_nt;
    }

    let mut known_spans: Vec<(String, usize)> = trunk_spans
        .iter()
        .filter(|(k, _)| !k.starts_with('?') && *k != "Unattributed")
        .map(|(k, &v)| (k.clone(), v))
        .collect();
    known_spans.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let backbone = if let Some((tr, _)) = known_spans.first() {
        tr.clone()
    } else {
        let mut all_spans: Vec<(String, usize)> =
            trunk_spans.iter().map(|(k, &v)| (k.clone(), v)).collect();
        all_spans.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        all_spans
            .first()
            .map(|x| x.0.clone())
            .unwrap_or_else(|| "Unknown".to_string())
    };

    let scope_suffix = match locus_scope {
        "pol" => "@pol",
        "env" => "@env",
        "partial" => "@partial",
        _ => "",
    };

    if trunk_spans.len() == 1 {
        let base = if locus_scope == "full" {
            format!("{}(Clonal)", backbone)
        } else {
            format!("{}(Clonal){}", backbone, scope_suffix)
        };
        return (base.clone(), backbone, base);
    }

    let mut donor_order: Vec<String> = Vec::new();
    let mut donor_genes: HashMap<String, Vec<String>> = HashMap::new();
    for s in segs {
        if s.trunk != backbone {
            if !donor_order.contains(&s.trunk) {
                donor_order.push(s.trunk.clone());
            }
            let entry = donor_genes.entry(s.trunk.clone()).or_default();
            if entry.last() != Some(&s.short_gene) {
                entry.push(s.short_gene.clone());
            }
        }
    }

    let donor_parts: Vec<String> = donor_order
        .into_iter()
        .map(|d| {
            let gs = &donor_genes[&d];
            format!("{}:{}", d, gs.join("+"))
        })
        .collect();

    let structural_formula = format!("{}[{}]", backbone, donor_parts.join(" | "));
    let descriptive_name = format!("{}{}", structural_formula, scope_suffix);
    (descriptive_name, backbone, structural_formula)
}

#[cfg(test)]
mod naming_tests {
    use super::*;

    fn seg(trunk: &str, start: usize, end: usize) -> SegmentCall {
        SegmentCall {
            seg_idx: 0,
            trunk: trunk.to_string(),
            hxb2_start: start,
            hxb2_end: end,
            span_nt: end - start + 1,
            ci_hxb2_start: start,
            ci_hxb2_end: end,
            short_gene: "pol".to_string(),
            identity_pct: 100.0,
        }
    }

    #[test]
    fn backbone_tiebreak_is_deterministic() {
        let a = seg("Unattributed", 800, 999);
        let b = seg("?E", 1000, 1199);
        for _ in 0..32 {
            let (_, bb1, _) = format_descriptive_name(&[a.clone(), b.clone()], "full");
            let (_, bb2, _) = format_descriptive_name(&[b.clone(), a.clone()], "full");
            assert_eq!(bb1, "?E");
            assert_eq!(bb2, "?E");
        }
    }
}

/// One-sided Fisher's exact test (`alternative = "greater"`) on a 2x2 contingency table `[[a, b], [c, d]]`.
pub fn fisher_exact_greater(a: usize, b: usize, c: usize, d: usize) -> f64 {
    let n = a + b + c + d;
    if n == 0 {
        return 1.0;
    }
    let r1 = a + b;
    let c1 = a + c;
    let max_a = r1.min(c1);

    let log_fact = |k: usize| -> f64 {
        if k <= 1 {
            0.0
        } else {
            // Stirling / Lanczos lgamma(k + 1)
            ln_gamma((k + 1) as f64)
        }
    };

    let log_denom = log_fact(n) - log_fact(r1) - log_fact(n - r1) - log_fact(c1) - log_fact(n - c1);
    let mut p_sum = 0.0f64;

    for x in a..=max_a {
        let y = r1 - x;
        let z = c1 - x;
        if n + x < r1 + c1 {
            continue;
        }
        let w = n + x - r1 - c1;
        let log_p = -log_denom - log_fact(x) - log_fact(y) - log_fact(z) - log_fact(w);
        p_sum += log_p.exp();
    }
    p_sum.clamp(0.0, 1.0)
}

/// Lanczos approximation for `ln(Gamma(z))`.
fn ln_gamma(z: f64) -> f64 {
    const COEFFS: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    let x = z - 1.0;
    let mut sum = COEFFS[0];
    for (i, &c) in COEFFS.iter().enumerate().skip(1) {
        sum += c / (x + i as f64);
    }
    let t = x + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + sum.ln()
}
