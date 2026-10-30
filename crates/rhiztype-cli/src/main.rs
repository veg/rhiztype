use clap::Parser;
use rayon::prelude::*;
use rhiztype_core::{RhizTypeEngine, ScreeningConfig, SequenceResult};
use std::fs;
use std::io::{self, Write};
use std::process::exit;
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

fn die(code: i32, msg: &str) -> ! {
    eprintln!("rhiztype: {}", msg);
    exit(code)
}

/// RFC-4180 quoted field; `guard` neutralises spreadsheet formula prefixes.
fn csv_quote(s: &str, guard: bool) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    if guard && s.starts_with(['=', '+', '-', '@']) {
        out.push('\'');
    }
    out.push_str(&s.replace('"', "\"\""));
    out.push('"');
    out
}

fn main() {
    let args = Args::parse();
    let fasta_path = match args.fasta.or(args.input_pos) {
        Some(p) => p,
        None => die(2, "please provide an input FASTA file via `<file.fasta>` or `--fasta <file.fasta>`"),
    };

    let engine_cfg = ScreeningConfig {
        gamma: args.gamma,
        min_span_nt: args.min_span,
        min_informative_snps: args.min_snps,
        p_fisher_thresh: args.p_thresh,
        mask_hypervariable_loops: !args.no_loop_mask,
    };
    if let Err(e) = engine_cfg.validate() {
        die(2, &format!("invalid configuration: {}", e));
    }

    let t0 = Instant::now();
    let bytes = fs::read(&fasta_path)
        .unwrap_or_else(|e| die(2, &format!("failed to read FASTA file {:?}: {}", fasta_path, e)));
    let raw_text = String::from_utf8(bytes)
        .unwrap_or_else(|_| die(2, &format!("FASTA file {:?} is not valid UTF-8", fasta_path)));

    let mut csv_file = args.out.as_ref().map(|p| {
        fs::File::create(p).unwrap_or_else(|e| die(2, &format!("cannot create CSV {:?}: {}", p, e)))
    });
    let mut json_file = args.json.as_ref().map(|p| {
        fs::File::create(p).unwrap_or_else(|e| die(2, &format!("cannot create JSON {:?}: {}", p, e)))
    });

    let mut records = RhizTypeEngine::parse_fasta(&raw_text);
    if records.is_empty() {
        die(1, "no FASTA records found");
    }
    let total_records = records.len();
    if args.limit > 0 && records.len() > args.limit {
        records.truncate(args.limit);
    }

    let engine = RhizTypeEngine::new();
    let cfg = engine_cfg;

    let results: Vec<SequenceResult> = records
        .par_iter()
        .map(|(hdr, seq)| engine.classify_sequence(hdr, seq, &cfg))
        .collect();

    let elapsed = t0.elapsed().as_secs_f64();

    if let Some(f) = csv_file.as_mut() {
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
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                csv_quote(&r.sequence_id, true),
                r.raw_length_nt,
                r.mapped_coding_nt,
                r.hxb2_span_start,
                r.hxb2_span_end,
                r.locus_scope,
                csv_quote(&r.alignment_mode, false),
                csv_quote(&r.descriptive_name, false),
                csv_quote(&r.structural_formula, false),
                r.backbone_trunk,
                r.n_segments,
                r.n_breakpoints,
                csv_quote(&r.hxb2_barcode, false),
                csv_quote(&r.ci_bracketed_barcode, false),
                csv_quote(&matches_joined, true)
            ));
        }
        let path = args.out.as_ref().unwrap();
        f.write_all((csv_lines.join("\n") + "\n").as_bytes())
            .unwrap_or_else(|e| die(2, &format!("failed to write CSV {:?}: {}", path, e)));
    }

    if let Some(f) = json_file.as_mut() {
        let json_str = serde_json::to_string_pretty(&results)
            .unwrap_or_else(|e| die(2, &format!("failed to serialize JSON: {}", e)));
        let path = args.json.as_ref().unwrap();
        f.write_all(json_str.as_bytes())
            .unwrap_or_else(|e| die(2, &format!("failed to write JSON {:?}: {}", path, e)));
    }

    let mut out = String::new();
    if !args.quiet {
        let n_label = if results.len() < total_records {
            format!("{} of {} records", results.len(), total_records)
        } else {
            format!("{} sequences", results.len())
        };
        out.push_str("================================================================================================\n");
        out.push_str(&format!(
            "RHIZTYPE HIV-1 GROUP M SUBTYPING & RECOMBINATION CLASSIFIER (N={}, {:.3}s)\n",
            n_label, elapsed
        ));
        out.push_str("================================================================================================\n");

        for r in &results {
            out.push_str(&format!("\n> {}\n", r.sequence_id));
            out.push_str(&format!(
                "  Alignment / Scope   : {} | {} (HXB2 {}..{}, {} mapped coding nt)\n",
                r.alignment_mode, r.locus_scope, r.hxb2_span_start, r.hxb2_span_end, r.mapped_coding_nt
            ));
            out.push_str(&format!("  Descriptive Formula : {}\n", r.descriptive_name));
            out.push_str(&format!("  HXB2 Barcode        : {}\n", r.hxb2_barcode));
            out.push_str(&format!("  CI-Bracketed Map    : {}\n", r.ci_bracketed_barcode));
            let match_strs: Vec<String> = r
                .matches
                .iter()
                .take(4)
                .map(|m| format!("{} [{}]", m.label, m.match_type))
                .collect();
            out.push_str(&format!("  Compendium Match    : {}\n", match_strs.join("; ")));
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
                out.push_str(&format!("  Verified Junctions  : {}\n", bp_strs.join("; ")));
            }
        }
    }
    if let Some(p) = &args.out {
        out.push_str(&format!("\n[+] Wrote CSV results ({} sequences) to {:?}\n", results.len(), p));
    }
    if let Some(p) = &args.json {
        out.push_str(&format!("[+] Wrote JSON report to {:?}\n", p));
    }

    // A closed pipe (e.g. `| head`) is not an error: files are already written.
    let stdout = io::stdout();
    let mut lock = stdout.lock();
    if let Err(e) = lock.write_all(out.as_bytes()).and_then(|_| lock.flush()) {
        if e.kind() != io::ErrorKind::BrokenPipe {
            die(2, &format!("failed to write to stdout: {}", e));
        }
    }
}
