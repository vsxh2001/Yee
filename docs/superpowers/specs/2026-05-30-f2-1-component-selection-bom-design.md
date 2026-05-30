# Filter Phase F2.1 — component selection + BOM — Design Spec

**ADR:** ADR-0112 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal

Second brick of the lumped-LC track (after F2.0 ideal LC ladder). Map each ideal
L/C value to a **real, purchasable standard component** (nearest IEC 60063
E-series value) and emit a **bill of materials**. Pure-data/math, WASM-safe, no
FDTD. Feeds the BOM panel in the UI and the tolerance analysis (F2.4, which sweeps
the realized-value response).

## Method

- **E-series** (IEC 60063): generate the standard preferred values for a decade
  and tile across decades. Support **E24** (±5 %) and **E96** (±1 %).
- **Nearest-value selection:** for an ideal value `x`, pick the E-series value
  minimizing `|log10(chosen) − log10(x)|` (log-nearest, the correct metric for a
  geometric series). Record the chosen value + the **deviation %**
  `(chosen − ideal)/ideal·100`. By construction the deviation is bounded by half
  the series ratio step (E24 ⇒ ≲ 5 %, E96 ⇒ ≲ 1 %).
- **BOM:** one line per ladder element (inductor or capacitor), with ideal value,
  chosen E-series value, deviation %, series, and a unit string (H / F shown in
  engineering units nH/pF). Group identical (kind, chosen-value) lines into a
  quantity. A `tolerance_pct` per line (the series tolerance) carries into F2.4.
- **Parasitics:** the BOM line carries OPTIONAL `esr_ohm`/`srf_hz` fields,
  defaulted `None` for the skeleton (a real vendor-parts library with measured
  ESR/SRF is a documented follow-on, F2.1b). Selection is value-only here.

## Changes (`crates/yee-filter/**` ONLY)

- New `crates/yee-filter/src/parts.rs`:
  - `pub enum ESeries { E24, E96 }` with `fn values_decade(&self) -> &[f64]` and
    `fn nearest(&self, x: f64) -> f64` (log-nearest across decades) +
    `fn tolerance_pct(&self) -> f64`.
  - `pub enum CompKind { Inductor, Capacitor }`
  - `pub struct BomLine { kind, ideal_value, chosen_value, deviation_pct,
    series: ESeries, tolerance_pct, qty, esr_ohm: Option<f64>, srf_hz: Option<f64> }`
  - `pub struct Bom { lines: Vec<BomLine> }` (+ serde) with a `fn total_parts()`.
  - `pub fn select_components(ladder: &LumpedLadder, series: ESeries) -> Bom`:
    for each `LcResonator`, add an inductor line (`l_henry`) and a capacitor line
    (`c_farad`), nearest-value-selected; group duplicates into `qty`.
- Re-export from `lib.rs`.

## DoD (machine-checkable; pure-math, NO FDTD)

1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-filter --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-filter` exit 0 — incl. a new gate `bom_001`:
   - **E-series correctness** (textbook anchors): `ESeries::E24.nearest(1.0)==1.0`,
     `nearest(4.5e-9)` ∈ {4.3e-9, 4.7e-9} picking the log-nearest, `nearest(1.05e3)`
     → 1.0e3 not 1.1e3 (log-nearest), and a handful more known cases; E96 picks a
     finer value. Every `nearest` result is an actual E-series member.
   - **Selection bound:** for the F2.0 cheb N=5 lumped ladder (E24), every chosen
     value is within the E24 quantization bound (≤ ~5.1 %) of its ideal.
   - **BOM completeness:** `total_parts() == 2·N` (an L + a C per resonator);
     grouping produces correct `qty` for the symmetric ladder (R1==R5, R2==R4 ⇒
     duplicate values merge).
   - Note: the gate validates the *selection*, NOT that the quantized response
     still passes the mask — that yield question is F2.4 (Monte-Carlo). Do NOT
     gate on quantized-response pass here.

## Out of scope

Vendor parts DB + measured ESR/SRF parasitics (F2.1b); the quantized-response /
yield analysis (F2.4); PCB footprints (F2.2); UI BOM panel; non-E24/E96 series.

## Why this next

Closed-form, pure-data, immediately validatable (E-series is deterministic), and
it produces the BOM + per-part tolerance the goal explicitly names — and the
tolerance metadata F2.4 needs. Depends only on F2.0's `LumpedLadder`.
