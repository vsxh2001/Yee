# App.2.2 — Low-pass stepped-impedance flow in the studio — Design Spec

**ADR:** ADR-0139 · **Date:** 2026-05-31 · **Status:** Accepted
**Vision:** `2026-05-31-ideal-filter-design-app-vision.md` §5 (response breadth — the
first **low-pass** capability in the visible app). Surfaces the F1.2.3
`dimension_stepped_impedance` engine (ADR-0137, shipped + Pozar-§8.6-gated but
un-surfaced) and makes the App.2.0 recommender's `SteppedImpedance` recommendation
routable to a real flow.

## Problem

The studio is **band-pass-only** (edge-coupled, hairpin, lumped — all band-pass). The
stepped-impedance **low-pass** dimensioner shipped (F1.2.3) but is un-surfaced; the
recommender recommends `SteppedImpedance` for low-pass but can only route to a
band-pass stand-in. The single biggest remaining capability gap is a **low-pass
response class** end-to-end in the app.

## Key insight (why this is integration, not new physics)

Both engine pieces exist and are validated: the low-pass magnitude response
(`yee_filter::lowpass_s21_squared`, the closed-form `|S21|²(Ω)` already used by the
band-pass `ideal_response` — for a low-pass filter evaluate at `Ω = f/f_c` with no
band-pass mapping) and the stepped-impedance dimensioner (`dimension_stepped_impedance`,
Pozar §8.6, ±0.02°). The **lumped flow is a parallel-response-path precedent** — a
`LumpedDesigned` + `design_lumped_from` + `lumped_*_stage`s + a LUMPED rail +
`StageCanvas` branching on `lumped_flow`. The low-pass stepped-Z flow mirrors that.

## Method

### 1. Engine API (`yee-filter`, small + strongly gated)

A public `ideal_response_lowpass(approx, order, cutoff_hz, freqs_hz) -> Vec<Complex64>`
that reuses the existing private `lowpass_s21_squared` at `Ω = f/f_c` (no band-pass
transform). This is the low-pass analogue of the existing public `ideal_response`.

### 2. Studio low-pass flow (`yee-studio-web`, mirror the lumped parallel flow)

- `Topology::SteppedImpedance` added; `Stage::rail(SteppedImpedance)` = a **stepped**
  rail (Spec / Technique / Synthesis / Layout / Export — the distributed five, no
  Components/Tolerance).
- A `SteppedLowpassDesigned` (mirror `LumpedDesigned`): from a **low-pass** spec
  (`Response::Lowpass`, `f0_hz` reused as the cutoff) it holds the prototype g-values,
  the stepped sections (`dimension_stepped_impedance`), the swept low-pass `|S21|`
  (`ideal_response_lowpass`) + the low-pass mask bands + PASS/FAIL, the board `Layout`,
  board size, and `dim_error`. `design_stepped_from(spec)` / `design_stepped()` +
  a `stepped` signal in `main.rs` recomputed on spec edit (mirror `lumped`).
- `stepped_synthesis_stage(stepped)`: the g-values + the swept low-pass `|S21|` vs the
  shaded low-pass mask + PASS/FAIL + the stepped-section table (Z, electrical length,
  width, length per section).
- `stepped_layout_stage(stepped)`: the stepped-Z board (the generic `Layout` SVG —
  reuse `board_svg`) + the section table + stackup.
- `StageCanvas`: `stepped_flow = topology() == SteppedImpedance` routes Synthesis /
  Layout to the stepped renderers (mirror `lumped_flow`).
- **Spec stage low-pass awareness:** when the active topology is `SteppedImpedance`,
  the Spec form labels the frequency as **Cutoff** and hides Fractional bandwidth
  (low-pass has no FBW); selecting the SteppedImpedance technique sets
  `spec.response = Response::Lowpass` (other techniques keep `Bandpass`).
- `technique_status`: `SteppedImpedance => Live(Topology::SteppedImpedance)`; card lit;
  `topology_label` / `topology_name` / `length_label` gain SteppedImpedance arms.
- Export: the stepped-Z Gerber / KiCad from the real `Layout` (reuse the generic export).

## Changes

- `crates/yee-filter/src/lib.rs` — `ideal_response_lowpass` (public) + a gate test.
- `crates/yee-studio-web/src/{engine.rs, stages.rs, main.rs}` — the parallel low-pass
  stepped-Z flow (mirror the lumped flow). `svg.rs` only if the board needs a tweak
  (it should not — `Layout` renders generically).

## DoD (machine-checkable)

1. **Engine gate (`yee-filter`, strong + non-vacuous):** `ideal_response_lowpass` —
   Butterworth: `|S21(f_c)|` = −3.01 dB ± 0.1 dB (the defining half-power cutoff);
   `|S21|` monotonically decreasing past `f_c`; deep stopband approaches the
   `−20·N·log10(f/f_c)` asymptote (assert a specific point, e.g. N=5 at `2 f_c` ≈
   −30 dB ± 1 dB). Chebyshev (ripple_db): in-band stays within the ripple bound and
   `|S21(f_c)|` = −ripple_db (equiripple edge). A constant response fails.
2. **Studio gate:** `dx build --platform web --release` EXIT 0; a NEW non-vacuous host
   test — `design_stepped_from(<low-pass demo spec>)` returns real stepped sections
   (from `dimension_stepped_impedance`, ≥ order sections, low-Z first) AND a swept
   `|S21|` that is ≈ −3 dB at the cutoff (Butterworth) — proving the SteppedImpedance
   card routes to the REAL low-pass engine, not a stub. Existing band-pass + lumped
   flows unregressed.
3. `cargo test -p yee-filter -p yee-studio-web` green; `cargo clippy ... -D warnings` +
   `cargo fmt --check` clean; `cargo check --workspace` green.

## Out of scope

Combline / interdigital; high-pass / band-stop; stepped-Z stopband-target auto-order;
EM verify (ADR-0133 wall untouched). Elliptic low-pass.

## Why

The first **low-pass** capability end-to-end in the visible app — the biggest remaining
response gap — surfacing two already-validated engines via the proven lumped
parallel-flow pattern. Strong low-pass-response gate + a non-vacuous routing gate.
