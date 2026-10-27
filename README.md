# RhizType — HIV-1 Group M Subtyping, Mosaic Deconvolution & Coordinate Nomenclature

**RhizType** is a high-throughput Rust command-line tool and client-side WebAssembly (WASM) browser application for screening **aligned or unaligned HIV-1 Group M FASTA sequences** (full-length genomes, partial `pol` PR-RT-IN surveillance amplicons, partial `env` gp120-gp41 sequences, or arbitrary coding fragments) for subtype assignment, mosaic breakpoint verification, and coordinate-resolved descriptive nomenclature.

- **Web Application (100% Client-Side WASM):** [https://veg.github.io/rhiztype/](https://veg.github.io/rhiztype/)
- **Source Repository:** [https://github.com/veg/rhiztype](https://github.com/veg/rhiztype)

---

## Key Capabilities

1. **Accepts Unaligned or Pre-Aligned FASTA Sequences Automatically**
   - **Pre-aligned inputs** mapped to HXB2 (`9,719 nt`), HXB2 coding span (`8,628 nt`, `HXB2 790–9417`), or the LANL 2023 HIV-1 Reference Compendium MSA (`11,531 columns`) are projected directly in $O(L)$ time.
   - **Unaligned sequences** (complete genomes or partial amplicons of any length $\ge 100\text{ nt}$) are automatically aligned to HXB2 coordinates using a two-stage **15-mer / 11-mer anchor seed hash + collinear sparse dynamic programming chain + banded Needleman-Wunsch** aligner against the Group M basal trunk ensemble.

2. **Scrubbed Basal Reference Panel (`13 Trunks + ?E`)**
   - Screens queries directly against the 13 scrubbed Group M basal reference trunks (`A1, A2, A6, B, C, D, F1, F2, G0, H, J, K, L`)—with the 12 historical reference-panel contaminations masked out—plus an **unpolarized ancestral root channel** and the unsampled `CRF01_AE` `vif..nef` (`?E`) trunk (`HXB2 5097–9417`).
   - Resolves deep-branching second-generation CRFs (`CRF01_AE`, `CRF02_AG`, `CRF06_cpx`, `CRF11_cpx`, `CRF63_02A6`, etc.) into their elementary basal building blocks without recursive label nesting.

3. **Hypervariable Loop Shield & Local $2 \times 2$ Fisher Exact Breakpoint Verification**
   - Automatically masks length-polymorphic hypervariable loops (`gp120` V1–V5 and `gag` p6 PTAP duplication hotspots) before segmentation so alignment indel noise never triggers false mosaic calls, while preserving all coding positions across `pol` (`PR`, `RT`, `IN`).
   - Evaluates every candidate recombination junction via a one-sided $2 \times 2$ Fisher's exact test ($p < 10^{-4}$ default) and a minimum informative-SNP support floor ($n_{11}, n_{22} \ge 5$), reporting exact informative-site plateau confidence intervals (`[CI:low-high]`).

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
# Screen example unaligned HIV-1 genomes (HXB2, CRF01_AE, CRF02_AG, CRF06_cpx)
./target/release/rhiztype docs/examples/demos.fasta

# Screen unaligned partial env sequences and export CSV + JSON reports
./target/release/rhiztype docs/examples/partial_env.fasta \
  --out env_results.csv \
  --json env_results.json

# Customize statistical stringency or disable DRM/loop masking
./target/release/rhiztype input.fasta \
  --p-fisher 1e-4 \
  --min-snps 5 \
  --threads 8 \
  --out results.csv
```

### Example Output

```text
================================================================================================
RHIZTYPE HIV-1 GROUP M SUBTYPING & RECOMBINATION CLASSIFIER (N=4 sequences, 0.005s)
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
  HXB2 Barcode        : A1:790-5096|?E:5097-9417
  CI-Bracketed Map    : A1:790-5096[CI:790-5151]|?E:5097-9417[CI:5097-9417]
  Compendium Match    : CRF01_AE [Exact Structural Formula]
  Verified Junctions  : HXB2 5097 (A1->?E, CI 5097-5151, p=3.89e-212)

> Ref.02_AG.NG.x.IBNG.L39106_unaligned
  Alignment / Scope   : Auto-aligned (15-mer anchor + banded NW to HXB2) | full (HXB2 790..9417, 8563 mapped coding nt)
  Descriptive Formula : A1[G0:gag..RT+RT..vif+gp41..nef]
  HXB2 Barcode        : A1:790-2237|G0:2238-3182|A1:3183-4119|G0:4120-5516|A1:5517-8159|G0:8160-9417
  CI-Bracketed Map    : A1:790-2237[CI:790-2241]|G0:2238-3182[CI:2237-3203]|A1:3183-4119[CI:3182-4128]|G0:4120-5516[CI:4119-5519]|A1:5517-8159[CI:5516-8227]|G0:8160-9417[CI:8159-9417]
  Compendium Match    : CRF02_AG [Homologous Mosaic (91.6% concordance, 4/5 breakpoints)]
```

---

## Rebuilding the WebAssembly (WASM) Bundle

To recompile the browser WebAssembly package into `docs/pkg/`:

```bash
rustup target add wasm32-unknown-unknown
wasm-pack build crates/rhiztype-wasm --target web --out-dir ../../docs/pkg --release
rm -f docs/pkg/.gitignore
```

To test the web application locally:

```bash
python3 -m http.server 8080 --directory docs
# Open http://localhost:8080 in your browser
```

---

## License

MIT License — Viral Evolution Group (VEG), Temple University.
