use clap::Parser;
use rayon::prelude::*;
use rhiztype_core::{RhizTypeEngine, ScreeningConfig, SequenceResult};
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(
    name = "rhiztype",
    version = "0.1.0",
    about = "RhizType: Fast HXB2-anchored HIV-1 Group M subtyping and mosaic deconvolution (aligned & unaligned FASTA)"
)]
struct Args {
    /// Input FASTA file (aligned or unaligned, full-length or partial pol/env)
    #[arg(short, long)]
    fasta: Option<PathBuf>,

    /// Positional FASTA path if --fasta is omitted
    #[arg(index = 1)]
    input_pos: Option<PathBuf>,

    /// Output CSV file path for tabular results
    #[arg(short, long)]
    out: Option<PathBuf>,

    /// Output JSON file path for full detailed segment/breakpoint report
    #[arg(long)]
    json: Option<PathBuf>,

    /// Maximum number of sequences to process (0 = all)
    #[arg(short, long, default_value_t = 0)]
    limit: usize,

    /// Template-switch penalty gamma (default 7.0 if omitted)
    #[arg(long)]
    gamma: Option<f64>,

    /// Minimum segment length in HXB2 nucleotides (auto-scales to 120 full / 150 partial if omitted)
    #[arg(long)]
    min_span: Option<usize>,

    /// Minimum subtype-informative SNPs per segment
    #[arg(long, default_value_t = 5)]
    min_snps: usize,

    /// Fisher exact test p-value threshold
    #[arg(long, default_value_t = 1e-4)]
    p_thresh: f64,

    /// Disable masking of hypervariable gp120 V1-V5 loops and gag p6 PTAP indel zone
    #[arg(long, default_value_t = false)]
    no_loop_mask: bool,

    /// Quiet mode (suppress per-sequence console report)
    #[arg(short, long, default_value_t = false)]
    quiet: bool,
}

fn main() {
    let args = Args::parse();
    let fasta_path = args
        .fasta
        .or(args.input_pos)
        .expect("Please provide an input FASTA file via `<file.fasta>` or `--fasta <file.fasta>`");

    let t0 = Instant::now();
    let raw_text = fs::read_to_string(&fasta_path)
        .unwrap_or_else(|e| panic!("Failed to read FASTA file {:?}: {}", fasta_path, e));

    let mut records = RhizTypeEngine::parse_fasta(&raw_text);
    if args.limit > 0 && records.len() > args.limit {
        records.truncate(args.limit);
    }

    let engine = RhizTypeEngine::new();
    let cfg = ScreeningConfig {
        gamma: args.gamma,
        min_span_nt: args.min_span,
        min_informative_snps: args.min_snps,
        p_fisher_thresh: args.p_thresh,
        mask_hypervariable_loops: !args.no_loop_mask,
    };

    let results: Vec<SequenceResult> = records
        .par_iter()
        .map(|(hdr, seq)| engine.classify_sequence(hdr, seq, &cfg))
        .collect();

    let elapsed = t0.elapsed().as_secs_f64();

    if !args.quiet {
        println!("================================================================================================");
        println!(
            "RHIZTYPE HIV-1 GROUP M SUBTYPING & RECOMBINATION CLASSIFIER (N={} sequences, {:.3}s)",
            results.len(),
            elapsed
        );
        println!("================================================================================================");

        for r in &results {
            println!("\n> {}", r.sequence_id);
            println!(
                "  Alignment / Scope   : {} | {} (HXB2 {}..{}, {} mapped coding nt)",
                r.alignment_mode, r.locus_scope, r.hxb2_span_start, r.hxb2_span_end, r.mapped_coding_nt
            );
            println!("  Descriptive Formula : {}", r.descriptive_name);
            println!("  HXB2 Barcode        : {}", r.hxb2_barcode);
            println!("  CI-Bracketed Map    : {}", r.ci_bracketed_barcode);
            let match_strs: Vec<String> = r
                .matches
                .iter()
                .take(4)
                .map(|m| format!("{} [{}]", m.label, m.match_type))
                .collect();
            println!("  Compendium Match    : {}", match_strs.join("; "));
            if !r.breakpoints.is_empty() {
                let bp_strs: Vec<String> = r
                    .breakpoints
                    .iter()
                    .map(|b| {
                        format!(
                            "HXB2 {} ({}->{}, CI {}-{}, p={:.2e})",
                            b.bp_hxb2, b.left_trunk, b.right_trunk, b.ci_hxb2_low, b.ci_hxb2_high, b.p_fisher
                        )
                    })
                    .collect();
                println!("  Verified Junctions  : {}", bp_strs.join("; "));
            }
        }
    }

    if let Some(out_csv) = args.out {
        let mut csv_lines = Vec::with_capacity(results.len() + 1);
        csv_lines.push(
            "sequence_id,raw_length_nt,mapped_coding_nt,hxb2_span_start,hxb2_span_end,locus_scope,alignment_mode,descriptive_name,structural_formula,backbone_trunk,n_segments,n_breakpoints,hxb2_barcode,ci_bracketed_barcode,compendium_matches".to_string(),
        );
        for r in &results {
            let matches_joined = r
                .matches
                .iter()
                .map(|m| format!("{} ({})", m.label, m.descriptive_name))
                .collect::<Vec<_>>()
                .join("; ");
            csv_lines.push(format!(
                "\"{}\",{},{},{},{},{},\"{}\",\"{}\",\"{}\",{},{},{},\"{}\",\"{}\",\"{}\"",
                r.sequence_id.replace('"', "'"),
                r.raw_length_nt,
                r.mapped_coding_nt,
                r.hxb2_span_start,
                r.hxb2_span_end,
                r.locus_scope,
                r.alignment_mode,
                r.descriptive_name,
                r.structural_formula,
                r.backbone_trunk,
                r.n_segments,
                r.n_breakpoints,
                r.hxb2_barcode,
                r.ci_bracketed_barcode,
                matches_joined.replace('"', "'")
            ));
        }
        fs::write(&out_csv, csv_lines.join("\n") + "\n")
            .unwrap_or_else(|e| panic!("Failed to write CSV {:?}: {}", out_csv, e));
        println!("\n[+] Wrote CSV results ({} sequences) to {:?}", results.len(), out_csv);
    }

    if let Some(out_json) = args.json {
        let json_str = serde_json::to_string_pretty(&results).expect("Failed to serialize JSON");
        fs::write(&out_json, json_str)
            .unwrap_or_else(|e| panic!("Failed to write JSON {:?}: {}", out_json, e));
        println!("[+] Wrote JSON report to {:?}", out_json);
    }
}
