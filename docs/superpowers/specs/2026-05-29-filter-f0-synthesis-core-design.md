# Filter Phase F0 — synthesis core (`yee-synth` + `yee-filter`) — Design Spec

**Phase:** F0 (filter roadmap)
**ADR:** ADR-0084
**Date:** 2026-05-29
**Status:** Accepted

---

## 1. Goal

The minimal end-to-end filter pipe: `FilterSpec → synthesized prototype +
coupling matrix → ideal response → spec-mask pass/fail`, exposed as
`yee filter synth`. Pure math; no EM; no new heavy dependency. Establishes the
data model every later phase plugs into.

## 2. `yee-synth` — synthesis math (no EM, no I/O)

### 2.1 Butterworth prototype (maximally flat)
`g0 = 1`; `g_k = 2·sin((2k−1)·π / (2N))`, k = 1..N; `g_{N+1} = 1`.

### 2.2 Chebyshev prototype (equi-ripple), ripple `L_Ar` dB, order N
(Pozar §8.3, eq 8.53):
```
β  = ln( coth( L_Ar / 17.37 ) )
γ  = sinh( β / (2N) )
a_k = sin( (2k−1)·π / (2N) ),           k = 1..N
b_k = γ² + sin²( k·π / N ),             k = 1..N
g0  = 1
g1  = 2·a_1 / γ
g_k = (4·a_{k−1}·a_k) / (b_{k−1}·g_{k−1}),  k = 2..N
g_{N+1} = 1                       (N odd)
g_{N+1} = coth²( β / 4 )          (N even)
```

### 2.3 Order estimation (optional helper)
Given ripple/RL + required rejection `A_s` dB at stopband ratio `Ω_s = ω_s/ω_c`:
Butterworth `N ≥ log10((10^{A_s/10}−1)/(10^{L_Ar/10}−1)) / (2·log10 Ω_s)`;
Chebyshev `N ≥ acosh(√((10^{A_s/10}−1)/(10^{L_Ar/10}−1))) / acosh(Ω_s)`.

### 2.4 Lowpass→bandpass transform
Centre `ω0 = √(ω1 ω2)`, fractional bandwidth `FBW = (ω2−ω1)/ω0`. Map prototype
Ω → `(1/FBW)·(ω/ω0 − ω0/ω)`.

### 2.5 Coupling coefficients + external Q (all-pole, synchronous)
```
k_{i,i+1} = FBW / √( g_i · g_{i+1} ),   i = 1..N−1
Qe_in  = g0·g1 / FBW
Qe_out = g_N·g_{N+1} / FBW
```
Normalized N×N coupling matrix `M` (synchronous → zero diagonal):
`M[i][i+1] = M[i+1][i] = 1/√(g_i g_{i+1})`; all other entries 0.

### 2.6 Public API (sketch)
```rust
pub enum Approximation { Butterworth, Chebyshev { ripple_db: f64 } }
pub struct Prototype { pub g: Vec<f64> }            // g[0]=g0 .. g[N+1]
pub fn prototype(approx: Approximation, order: usize) -> Prototype;
pub fn min_order(approx: Approximation, rejection_db: f64, omega_s: f64) -> usize;
pub struct CouplingDesign {                          // §2.5 outputs
    pub k: Vec<f64>, pub qe_in: f64, pub qe_out: f64, pub m: Vec<Vec<f64>>,
}
pub fn coupling_design(proto: &Prototype, fbw: f64) -> CouplingDesign;
```
`#![forbid(unsafe_code)]`, `#![warn(missing_docs)]`. Dep: `yee-core`, `nalgebra`.

## 3. `yee-filter` — data model + ideal response + flow scaffold

```rust
pub enum Response { Lowpass, Highpass, Bandpass, Bandstop }
pub struct FilterSpec {              // serde
    pub response: Response,
    pub approximation: Approximation, // re-exported from yee-synth
    pub f0_hz: f64, pub fbw: f64,     // or f1/f2; bandpass
    pub order: Option<usize>,         // None → min_order from mask
    pub z0_ohm: f64,
    pub mask: SpecMask,
}
pub struct SpecMask {                // passband + stopband points
    pub passband_ripple_db: f64, pub return_loss_db: f64,
    pub stopband: Vec<(f64 /*hz*/, f64 /*min reject dB*/)>,
}
pub struct CouplingMatrix { pub m: Vec<Vec<f64>>, pub qe_in: f64, pub qe_out: f64 }
pub enum Topology { CoupledResonator, /* future: Ladder, Iris, ... */ }
pub struct FilterProject {           // the persisted design document (serde)
    pub spec: FilterSpec, pub prototype: Prototype,
    pub coupling: CouplingMatrix, pub topology: Topology,
}
pub fn synthesize(spec: &FilterSpec) -> FilterProject;          // Stages 1–2
pub fn ideal_response(proj: &FilterProject, freqs_hz: &[f64]) -> Vec<Complex64>; // S21
pub fn check_mask(proj: &FilterProject, freqs_hz: &[f64]) -> MaskReport;         // Stage-3 gate
```
Ideal bandpass `|S21|²(ω)` from the **closed-form** transfer function applied to
the bandpass-mapped Ω (§2.4): Chebyshev `1/(1+ε²T_N²(Ω))`, `ε=√(10^{L_Ar/10}−1)`,
`T_N`=Chebyshev poly; Butterworth `1/(1+Ω^{2N})`. `S11` from `|S11|²=1−|S21|²`
(lossless). Dep: `yee-core`, `yee-synth`, `yee-io` (Touchstone), `serde`.

## 4. `yee-cli` — `yee filter synth <spec.toml>`
Parse a `FilterSpec` TOML; `synthesize`; print prototype g-values, coupling
matrix, Qe; sweep `ideal_response` over a band; write Touchstone (`yee-io`);
print the `check_mask` verdict; exit 0 on pass, 1 on mask fail. A new `Filter`
subcommand with a `Synth { spec, output }` variant. Add an example
`crates/yee-cli/examples` or a fixture spec under the crate.

## 5. Definition of Done (machine-checkable)

1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-synth -p yee-filter -p yee-cli --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-synth -p yee-filter -p yee-cli` exit 0 (all fast; no EM).
4. **`synth-001`** (yee-synth test): computed g-values match published tables to
   `≤ 1e-3` absolute:
   - Butterworth N=3 → `[1.0, 2.0, 1.0]`; N=5 → `[0.6180, 1.6180, 2.0000, 1.6180, 0.6180]`.
   - Chebyshev 0.5 dB N=3 → `g1,g2,g3 = 1.5963, 1.0967, 1.5963`, `g4=1.0`.
   - Chebyshev 0.5 dB N=5 → `1.7058, 1.2296, 2.5408, 1.2296, 1.7058`, `g6=1.0`.
   - Chebyshev 3.0 dB N=3 → `3.3487, 0.7117, 3.3487`, `g4=1.0`.
   - Chebyshev 0.5 dB N=4 → `g5 = coth²(β/4)` ≈ `1.9841` (even-order load check).
5. **`synth-002`** (yee-synth test): for a worked example (e.g. Chebyshev 0.5 dB
   N=3, FBW=0.10) the coupling coefficients `k_{12}=k_{23}` and `Qe` match the
   §2.5 closed form recomputed independently in the test (and, where a published
   Hong-Lancaster value exists, within `≤ 1e-3`).
6. **`filt-001`** (yee-filter test): `synthesize` a Chebyshev 0.5 dB bandpass,
   then `check_mask` over a swept band returns PASS — passband ripple ≤ 0.5 dB,
   in-band RL ≥ the spec, and rejection ≥ mask at the stopband points. A
   deliberately-too-low order returns FAIL (negative control).
7. `yee filter synth <fixture-spec.toml>` exits 0, prints the matrix, writes a
   readable Touchstone, and reports PASS for a satisfiable spec.

The `synth-001`/`synth-002`/`filt-001` gates live as `#[test]`s in
`crates/yee-synth/tests/` and `crates/yee-filter/tests/` — this satisfies the
CLAUDE.md §4 contract ("a published-benchmark validation case in
`crates/<crate>/tests/`"). Registering them in the `yee-validation` aggregator
(so they appear in `yee validate --list`) is a small follow-on, **Phase F0.1**,
kept out of F0 to hold the lane to the two new crates + the CLI and avoid a new
`yee-validation → yee-synth/yee-filter` dependency edge in the skeleton.

## 6. Out of scope (later phases)
Coupling-matrix→S realization; elliptic/Cameron synthesis; EM; layout; export;
GUI; `yee-validation` aggregator registration (Phase F0.1). F0 is math + data
model + CLI only.
