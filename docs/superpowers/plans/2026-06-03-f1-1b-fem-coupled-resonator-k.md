# F1.1b coupling-coefficient extraction via FEM driven-sweep — implementation plan

**Spec:** [2026-06-03-f1-1b-fem-coupled-resonator-k-design.md](../specs/2026-06-03-f1-1b-fem-coupled-resonator-k-design.md)
**ADR:** [ADR-0155](../../src/decisions/0155-f1-1b-fem-coupled-resonator-k.md)
**Seed:** `spike/fem-coupled-k-probe` (`933940f`) — `crates/yee-fem/tests/coupled_k_probe.rs` (GO probe).

## Brick K1 — coupled-resonator-k API + `fem-coupling-001` gate

1. Confirm dep direction: `yee-filter` must NOT depend on `yee-fem` (so `yee-fem` may dep
   `yee-filter` for `extract_coupling`). If acyclic, add `yee-filter` to `yee-fem` deps (real dep if
   the API lives in `src` and calls `extract_coupling`; else dev-dep + inline the peak-finder).
2. New `crates/yee-fem/src/coupled_resonator_k.rs` (promote the probe's geometry+sweep helpers from
   `coupled_k_probe.rs`): `pub struct CoupledResonatorGeom { trace_w, gap_s, sub_h, eps_r, f0_hz, box_w, box_h }`
   + `pub struct CoupledKResult { f_lo_hz, f_hi_hz, k_fem, k_imp_ref, k_eps_ref, peaks_resolvable, valley_db }`
   + `pub fn coupled_resonator_k(geom: &CoupledResonatorGeom, n_pts: usize) -> Result<CoupledKResult, Error>`.
   Build the two-λ_g/2-resonator mesh via `layered_microstrip_filter_mesh` + `TraceRect`; two weak
   gap-coupled feeds via `microstrip_port_numerical_at`; `sweep_matrix` + `with_coupled_whitney(true)`;
   `yee_filter::extract_coupling` on |S21|(f). Re-export from `lib.rs`.
3. Fast unit test (non-ignored, debug-safe): geometry well-formed; `k_imp` via
   `coupling_coefficient(&coupled_microstrip(...))` and `k_eps` finite + positive; the two-k
   references agree to ~order (the probe's premise). No FEM solve.
4. Gate `crates/yee-fem/tests/fem_coupling_001.rs` (`#[ignore]`'d + `--release`): call
   `coupled_resonator_k`; assert `peaks_resolvable` + `valley_db ≤ shallower_peak_db − 6.0` (a real
   margin) + `|k_fem − k_imp_ref| / k_imp_ref ≤ 0.30`; print `k_fem`, `k_imp_ref`, `k_eps_ref`, the
   peak freqs, and the full |S21|(f). Keep the probe's geometry (W=1, S=2 mm, h=1, ε_r=4.4, f0=2.4 GHz).
5. CI: add a `fem-coupling-001` step to the existing `fem-eigen` `--release` gate job in
   `.github/workflows/ci.yml` (`libfontconfig1-dev` already there).
6. Verify (boxed): `cargo fmt --check`; `cargo clippy -p yee-fem --all-targets -- -D warnings`;
   `cargo test -p yee-fem --lib` (fast unit); `cargo test -p yee-fem --release --test fem_coupling_001
   -- --ignored --nocapture` → exit 0, prints k_fem≈0.048, two peaks. Boxed `scripts/yee-box.sh`.

## Brick K2 — k-vs-gap monotonicity (after K1 merge)

7. A `fem-coupling-002` gate: sweep S ∈ {1.5, 2.0, 3.0} mm, assert `k_fem(S)` monotonic-decreasing +
   each within 30 % of `coupling_coefficient(S)`. Heavier (3 sweeps) → release gate.

## Brick K3 — Qe (deferred, later increment)

## Dispatch

- **K1 = one agent**, worktree off `931...`/main, base the probe seed `933940f` (it has the working
  geometry + sweep), lane `crates/yee-fem/src/**`, `crates/yee-fem/tests/{coupled_k_probe.rs →
  promote, fem_coupling_001.rs}`, `crates/yee-fem/Cargo.toml`, `.github/workflows/ci.yml`.
- Reviewer (never self-review) after K1; I verify boxed + merge `--no-ff`.
- K2 after K1 merges.

## CI / discipline

- All heavy gates `#[ignore]`'d + `--release` job; never the debug workspace test.
- Honesty: the two-peaks tripwire + k-tolerance are real measured quantities; do NOT weaken the
  tolerance. A NO-GO geometry prints its measurement (the probe is non-failing by design; K1's gate
  is the real PASS/FAIL).
