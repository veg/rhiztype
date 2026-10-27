use rhiztype_core::{RhizTypeEngine, ScreeningConfig};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct WasmRhizType {
    engine: RhizTypeEngine,
}

#[wasm_bindgen]
impl WasmRhizType {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            engine: RhizTypeEngine::new(),
        }
    }

    /// Screen a FASTA string (aligned or unaligned, full-length or partial pol/env)
    /// and return a JSON array of `SequenceResult` objects.
    pub fn screen_fasta(
        &self,
        fasta_text: &str,
        gamma: f64,
        min_span_nt: usize,
        min_informative_snps: usize,
        p_fisher_thresh: f64,
        mask_hypervariable_loops: bool,
    ) -> String {
        let records = RhizTypeEngine::parse_fasta(fasta_text);
        let cfg = ScreeningConfig {
            gamma: if gamma > 0.0 { Some(gamma) } else { None },
            min_span_nt: if min_span_nt > 0 {
                Some(min_span_nt)
            } else {
                None
            },
            min_informative_snps,
            p_fisher_thresh,
            mask_hypervariable_loops,
        };

        let results: Vec<_> = records
            .iter()
            .map(|(hdr, seq)| self.engine.classify_sequence(hdr, seq, &cfg))
            .collect();

        serde_json::to_string(&results).unwrap_or_else(|_| "[]".to_string())
    }

    /// Return the embedded 200-lineage reference compendium as a JSON array.
    pub fn get_compendium_json(&self) -> String {
        serde_json::to_string(&self.engine.compendium).unwrap_or_else(|_| "[]".to_string())
    }
}
