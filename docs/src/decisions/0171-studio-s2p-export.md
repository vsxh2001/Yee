# ADR-0171: Studio `.s2p` Touchstone export — designed-filter S-parameter portability (T8)

**Status:** Accepted
**Date:** 2026-06-06
**Related:** ADR-0130 (Dioxus `yee-studio-web` is THE studio), ADR-0163 (studio finite-Q response), ADR-0161
(CLI `write_s2p` + the lossless-S11-quadrature passivity fix), `yee_io::touchstone`,
[[project-filter-design-final-goal]] ("see the parameters at the end" / manufacturing + interchange out).

---

## Context

The studio's Export stage emits Gerber, KiCad, BOM/CPL (lumped) and a parameter sheet — but NOT the designed
filter's **S-parameters as a Touchstone `.s2p`**. That is the standard RF interchange format: a user can import
it into scikit-rf / ADS / AWR / a VNA comparison. The CLI already writes `.s2p` (`yee_cli::filter::write_s2p`
→ `yee_io::touchstone::write`); the studio has no equivalent, so a design can't leave the app as S-parameters.
This is a concrete portability gap and a commonly-expected filter-tool feature.

The pieces exist: the complex S-parameter models are already in the WASM-safe path — `ladder_s_params_lossy`
returns `(S11, S21)` for the lumped finite-Q response, and `ideal_response` returns the complex `S21` for the
distributed/coupling-matrix flow (with `S11` placed in lossless quadrature, the ADR-0161 passivity-correct
form). The Touchstone *renderer* also exists — `yee_io::touchstone::render(&File) -> Result<String>` — but it
is **private**; `write(path, file)` is just `render` + `std::fs::write`.

## Decision

1. **`yee-io`: expose the WASM-safe Touchstone string renderer.** Make `render` public as
   `pub fn to_string(file: &File) -> Result<String>` (keep `render` as its internal name or rename); `write`
   delegates to it (no format duplication — single source of truth for the Touchstone text). Pure string
   formatting, no `fs`/`Path` — WASM-safe.
2. **`yee-studio-web`: a `.s2p` download** in the lumped + distributed Export panels. An engine helper
   `s2p_string(z0, freqs, &[(Complex64, Complex64)], comment) -> Result<String, _>` builds a 2-port
   `touchstone::File` (`n_ports = 2`, `Format::RealImag`, `FreqUnit::Hz`, row-major `[S11, S21, S21, S11]`
   reciprocal+symmetric — mirroring the CLI's `write_s2p`) and calls `yee_io::touchstone::to_string`. The
   download buttons compute the `(S11, S21)` sweep:
   - **Lumped:** the **finite-Q realistic** response over the sweep grid via `ladder_s_params_lossy(&ladder, f,
     q_unloaded)` (the curve the user sees + a VNA would measure). File named `filter-lumped.s2p`.
   - **Distributed:** `ideal_response(&project, &freqs)` for `S21` + the **lossless `S11` quadrature**
     (`S11 = j·√(1−|S21|²)`, the ADR-0161 passivity-correct pair) so `to_string`'s passivity check accepts it.
     Named `filter.s2p`.
   The `.s2p` matches the displayed response (same model/grid), so it is coherent with the plot regardless of
   the JLCPCB export's auto-routed topology (the top-C export-coherence is the separate deferred T7 arc).

**Gate (`yee-studio-web` engine, fast, non-`#[ignore]`'d, NON-circular):** build the `(S11, S21)` sweep for a
lumped fixture, render via `s2p_string`, and assert (a) it **round-trips** through a Touchstone parse
(`yee_io::touchstone` parse on the string, or re-validate header `# Hz S RI R <z0>` + `n_ports = 2` + one data
row per frequency), and (b) the rendered values match the input `(S11, S21)` within float tolerance. Mirrors the
CLI's `.s2p` contract at the studio-engine layer.

## Consequences

- The studio gains S-parameter portability — a design leaves the app as a standard `.s2p`, importable into the
  RF-tool ecosystem. Studio↔CLI parity on `.s2p`.
- Scope T8: `crates/yee-io/src/touchstone.rs` (renderer visibility only) + `crates/yee-studio-web/src/{engine.rs,
  stages.rs}` + the studio gate. WASM-safe (the renderer is pure text; `yee-io` compiles to wasm32 — verify no
  native-only transitive; keep the `wasm-build` job green).
- **Not in scope:** an `.s1p`/multi-port path (2-port filter); routing the `.s2p` to the auto-routed top-C
  response (the T7 coherence arc — the `.s2p` matches the *displayed* response, which is the right contract);
  CLI `write_s2p` dedup onto the shared renderer (it already calls `write` = `render` + fs; an optional later tidy).

## References
- Renderer: `yee_io::touchstone::{render (→ pub to_string), write, File, Format, FreqUnit}`.
- Models: `yee_filter::{ladder_s_params_lossy, ideal_response}`; the lossless `S11` quadrature
  (`lossless_s_pair`, ADR-0161 / `yee_cli::filter`).
- Studio: `crates/yee-studio-web/src/{engine.rs (design_lumped_from, Designed/LumpedDesigned), stages.rs
  (export_lumped, export_distributed, download_btn/download_file)}`.
