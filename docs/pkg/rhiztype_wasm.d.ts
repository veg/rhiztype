/* tslint:disable */
/* eslint-disable */

export class WasmRhizType {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Return the embedded 200-lineage reference compendium as a JSON array.
     */
    get_compendium_json(): string;
    constructor();
    /**
     * Screen a FASTA string (aligned or unaligned, full-length or partial pol/env)
     * and return a JSON array of `SequenceResult` objects.
     */
    screen_fasta(fasta_text: string, gamma: number, min_span_nt: number, min_informative_snps: number, p_fisher_thresh: number, mask_hypervariable_loops: boolean, mask_pol_drm: boolean): string;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_wasmrhiztype_free: (a: number, b: number) => void;
    readonly wasmrhiztype_get_compendium_json: (a: number) => [number, number];
    readonly wasmrhiztype_new: () => number;
    readonly wasmrhiztype_screen_fasta: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number) => [number, number];
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
