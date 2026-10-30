# RhizType — HIV-1 Group M Subtyping, Mosaic Deconvolution & Coordinate Nomenclature

**RhizType** is a high-throughput Rust command-line tool and client-side WebAssembly (WASM) browser application for screening **aligned or unaligned HIV-1 Group M FASTA sequences** (full-length genomes, partial `pol` PR-RT-IN surveillance amplicons, partial `env` gp120-gp41 sequences, or arbitrary coding fragments) for subtype assignment, mosaic breakpoint verification, and coordinate-resolved descriptive nomenclature.

- **Web Application (100% Client-Side WASM):** [https://veg.github.io/rhiztype/](https://veg.github.io/rhiztype/)
- **Source Repository:** [https://github.com/veg/rhiztype](https://github.com/veg/rhiztype)

---

## Key Capabilities

1. **Accepts Unaligned or Pre-Aligned FASTA Sequences Automatically**
   - **Pre-aligned inputs** mapped to HXB2 (`9,719 nt`), HXB2 coding span (`8,628 nt`, `HXB2 790–9417`), or the LANL 2023 HIV-1 Reference Compendium MSA (`11,531 columns`) are projected directly in $O(L)$ time.
   - **Unaligned sequences** (complete genomes or partial amplicons; at least $100\text{ nt}$ of usable coding positions are required, counted as unambiguous mapped coding positions before hypervariable-loop masking, not raw input length) are automatically aligned to HXB2 coordinates using a two-stage **15-mer / 11-mer anchor seed hash + collinear sparse dynamic programming chain + banded Needleman-Wunsch** aligner against the Group M basal trunk ensemble.

2. **Scrubbed Basal Reference Panel (`13 Trunks + ?E`)**
   - Screens queries directly against the 13 scrubbed Group M basal reference trunks (`A1, A2, A6, B, C, D, F1, F2, G0, H, J, K, L`)—with the 12 historical reference-panel contaminations masked out—plus an **unpolarized ancestral root channel** and the unsampled `CRF01_AE` `vif..nef` (`?E`) trunk (`HXB2 5097–9417`).
   - Resolves deep-branching second-generation CRFs (`CRF01_AE`, `CRF02_AG`, `CRF06_cpx`, `CRF11_cpx`, `CRF63_02A6`, etc.) into their elementary basal building blocks without recursive label nesting.

3. **Hypervariable Loop Shield & Local $2 \times 2$ Fisher Exact Breakpoint Verification**
   - Automatically masks length-polymorphic hypervariable loops (`gp120` V1–V5 and `gag` p6 PTAP duplication hotspots) before segmentation so alignment indel noise never triggers false mosaic calls, while preserving all coding positions across `pol` (`PR`, `RT`, `IN`).
   - Evaluates every candidate recombination junction via a one-sided $2 \times 2$ Fisher's exact test ($p < 10^{-4}$ default) and a minimum informative-SNP support floor ($n_{11}, n_{22} \ge 5$), reporting informative-site plateau brackets (`[CI:low-high]`: nearest flanking informative sites; not a coverage-probability CI; may be one-sided or zero-width at junctions with `?E`/unattributed segments). With perfect separation, the default $p < 10^{-4}$ requires roughly $\ge 8$ informative sites per side (5|5 gives $p \approx 4\times10^{-3}$), so the $n_{11}, n_{22} \ge 5$ floor is necessary but not sufficient. Junctions involving the unattributed (ghost) state test novel-mutation density, with the known trunk's matches as the comparison row, rather than the pairwise informative-SNP count.

4. **Dual-Layer Coordinate Nomenclature & 200-Lineage Compendium Matching**
   - Outputs both human-readable **Descriptive Structural Formulas** (e.g., `A1[?E:vif..nef]`, `C[B:gag.CA..RT+vif..vpr]`, `B[F1:RT]@pol`) and machine-readable **HXB2 Coordinate Barcodes** (`A1:790-5096|?E:5097-9417`), cross-referenced against 200 cataloged lineages (162 LANL subtypes/CRFs + 38 novel epidemiologically circulating full-genome, `pol`, and `env` lineages).

---

## Workspace Architecture

```text
rhiztype/
├── Cargo.toml                  # Rust workspace manifest
├── crates/
│   ├── rhiztype-core/          # Core alignment, Viterbi PA-FF segmentation, Fisher exact & compendium engine
│   ├── rhiztype-cli/           # Parallel (Rayon) native command-line binary (`rhiztype`)
│   └── rhiztype-wasm/          # WebAssembly (`wasm-bindgen`) bindings for browser execution
└── docs/                       # Static GitHub Pages web app + compiled WASM bundle + example FASTAs
    ├── index.html
    ├── pkg/
    └── examples/
```

---

## CLI Installation & Usage

### Build from Source

Requires Rust (`1.75+`):

```bash
git clone https://github.com/veg/rhiztype.git
cd rhiztype
cargo build --release -p rhiztype-cli
```

The optimized binary is placed at `./target/release/rhiztype` (and can be installed into `~/.cargo/bin` via `cargo install --path crates/rhiztype-cli`).

### Screen Sequences (Unaligned or Aligned FASTA)

```bash
# Screen example HIV-1 genomes (5 isolates: HXB2, CRF01_AE, CRF02_AG, CRF06_cpx, CRF07_BC 98CN009)
./target/release/rhiztype docs/examples/demos.fasta

# Screen unaligned partial env sequences and export CSV + JSON reports
./target/release/rhiztype docs/examples/partial_env.fasta \
  --out env_results.csv \
  --json env_results.json

# Customize statistical stringency or disable loop masking
./target/release/rhiztype input.fasta \
  --p-thresh 1e-4 \
  --min-snps 5 \
  --out results.csv

# Set parallelism with the Rayon environment variable
RAYON_NUM_THREADS=8 ./target/release/rhiztype input.fasta --out results.csv
```

| Flag | Meaning |
|---|---|
| `--fasta`/`-f` or positional | Input FASTA (aligned or unaligned) |
| `--out`/`-o` | CSV report path |
| `--json` | JSON report path |
| `--limit`/`-l` | Maximum number of records to process |
| `--gamma` | Viterbi switch penalty (default 7.0) |
| `--min-span` | Minimum segment span in nt (auto: 120 full / 150 partial) |
| `--min-snps` | Minimum informative SNPs per segment (5) |
| `--p-thresh` | Fisher exact p-value threshold (1e-4) |
| `--no-loop-mask` | Disable hypervariable loop masking |
| `--quiet`/`-q` | Suppress console report |

### Example Output

```text
================================================================================================
RHIZTYPE HIV-1 GROUP M SUBTYPING & RECOMBINATION CLASSIFIER (N=5 sequences, 0.224s)
================================================================================================

> Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455_unaligned
  Alignment / Scope   : Pre-aligned (HXB2 9,719 nt) | full (HXB2 790..9417, 8628 mapped coding nt)
  Descriptive Formula : B(Clonal)
  HXB2 Barcode        : B:790-9417
  CI-Bracketed Map    : B:790-9417[CI:790-9417]
  Compendium Match    : Subtype_B [Exact Structural Formula]

> Ref.01_AE.TH.90.CM240.U54771_unaligned
  Alignment / Scope   : Auto-aligned (15-mer anchor + banded NW to HXB2) | full (HXB2 790..9417, 8567 mapped coding nt)
  Descriptive Formula : A1[?E:vif..nef]
  HXB2 Barcode        : A1:790-5150|?E:5151-9417
  CI-Bracketed Map    : A1:790-5150[CI:790-5151]|?E:5151-9417[CI:5151-9417]
  Compendium Match    : CRF01_AE [Exact Structural Formula]; CRF100_01C [Exact Structural Formula]; CRF114_0155 [Exact Structural Formula]; CRF101_01B [Homologous Mosaic (95.1% concordance, 1/2 breakpoints)]
  Verified Junctions  : HXB2 5151 (A1->?E, CI 5151-5151, p=3.89e-216)

> Ref.02_AG.NG.x.IBNG.L39106_unaligned
  Alignment / Scope   : Auto-aligned (15-mer anchor + banded NW to HXB2) | full (HXB2 790..9417, 8563 mapped coding nt)
  Descriptive Formula : A1[G0:gag..RT+RT..vif+tat1..gp120+nef]
  HXB2 Barcode        : A1:790-2237|G0:2238-3194|A1:3195-4119|G0:4120-5516|A1:5517-5989|G0:5990-6332|A1:6333-8909|G0:8910-9417
  CI-Bracketed Map    : A1:790-2237[CI:790-2241]|G0:2238-3194[CI:2237-3203]|A1:3195-4119[CI:3194-4128]|G0:4120-5516[CI:4119-5519]|A1:5517-5989[CI:5516-5992]|G0:5990-6332[CI:5989-6341]|A1:6333-8909[CI:6332-8912]|G0:8910-9417[CI:8909-9417]
  Compendium Match    : CRF02_AG [Homologous Mosaic (83.8% concordance, 4/5 breakpoints)]
  Verified Junctions  : HXB2 2238 (A1->G0, CI 2237-2241, p=3.42e-12); HXB2 3195 (G0->A1, CI 3194-3203, p=6.04e-12); HXB2 4120 (A1->G0, CI 4119-4128, p=1.49e-11); HXB2 5517 (G0->A1, CI 5516-5519, p=4.46e-7); HXB2 5990 (A1->G0, CI 5989-5992, p=9.39e-6); HXB2 6333 (G0->A1, CI 6332-6341, p=7.22e-8); HXB2 8910 (A1->G0, CI 8909-8912, p=3.21e-11)

> Ref.06_cpx.AU.96.BFP90.AF064699_unaligned
  Alignment / Scope   : Auto-aligned (15-mer anchor + banded NW to HXB2) | full (HXB2 790..9417, 8576 mapped coding nt)
  Descriptive Formula : G0[A1:gag.MA..gag.CA+RT..IN | K:RT | J:IN..tat1+gp41..nef]
  HXB2 Barcode        : A1:790-1470|G0:1471-3209|K:3210-4136|A1:4137-4861|J:4862-6020|G0:6021-8261|J:8262-9417
  CI-Bracketed Map    : A1:790-1470[CI:790-1473]|G0:1471-3209[CI:1470-3224]|K:3210-4136[CI:3209-4138]|A1:4137-4861[CI:4136-4925]|J:4862-6020[CI:4861-6021]|G0:6021-8261[CI:6020-8895]|J:8262-9417[CI:8261-9417]
  Compendium Match    : CRF06_cpx [Homologous Mosaic (90.2% concordance, 5/6 breakpoints)]
  Verified Junctions  : HXB2 1471 (A1->G0, CI 1470-1473, p=5.46e-15); HXB2 3210 (G0->K, CI 3209-3224, p=1.23e-9); HXB2 4137 (K->A1, CI 4136-4138, p=6.65e-9); HXB2 4862 (A1->J, CI 4861-4925, p=3.92e-9); HXB2 6021 (J->G0, CI 6020-6021, p=9.29e-23); HXB2 8262 (G0->J, CI 8261-8895, p=3.84e-44)

> Ref.07_BC.CN.98.98CN009.AF286230_unaligned
  Alignment / Scope   : Auto-aligned (15-mer anchor + banded NW to HXB2) | full (HXB2 790..9417, 8579 mapped coding nt)
  Descriptive Formula : C[B:gag..PR+RT+vpr..gp120+nef]
  HXB2 Barcode        : C:790-1980|B:1981-2541|C:2542-2939|B:2940-3179|C:3180-5679|B:5680-6428|C:6429-8842|B:8843-9057|C:9058-9417
  CI-Bracketed Map    : C:790-1980[CI:790-1990]|B:1981-2541[CI:1980-2552]|C:2542-2939[CI:2541-2953]|B:2940-3179[CI:2939-3206]|C:3180-5679[CI:3179-5680]|B:5680-6428[CI:5679-6434]|C:6429-8842[CI:6428-8847]|B:8843-9057[CI:8842-9061]|C:9058-9417[CI:9057-9417]
  Compendium Match    : CRF07_BC [Homologous Mosaic (95.0% concordance, 6/6 breakpoints)]
  Verified Junctions  : HXB2 1981 (C->B, CI 1980-1990, p=2.58e-13); HXB2 2542 (B->C, CI 2541-2552, p=1.33e-12); HXB2 2940 (C->B, CI 2939-2953, p=3.71e-11); HXB2 3180 (B->C, CI 3179-3206, p=1.08e-14); HXB2 5680 (C->B, CI 5679-5680, p=4.38e-22); HXB2 6429 (B->C, CI 6428-6434, p=1.20e-25); HXB2 8843 (C->B, CI 8842-8847, p=1.53e-13); HXB2 9058 (B->C, CI 9057-9061, p=1.09e-11)
```

---

## Rebuilding the WebAssembly (WASM) Bundle

To recompile the browser WebAssembly package into `docs/pkg/`:

```bash
rustup target add wasm32-unknown-unknown
wasm-pack build crates/rhiztype-wasm --target web --out-dir ../../docs/pkg --release && rm -f docs/pkg/.gitignore
```

`docs/pkg/` is a committed binary artifact: rebuild and commit it whenever `crates/rhiztype-core` changes, since CI only uploads `docs/` and does not deploy a rebuilt bundle. Run `cargo test --workspace` before rebuilding.

To test the web application locally:

```bash
python3 -m http.server 8080 --directory docs
# Open http://localhost:8080 in your browser
```

---

## License

MIT License — Viral Evolution Group (VEG), Temple University.
