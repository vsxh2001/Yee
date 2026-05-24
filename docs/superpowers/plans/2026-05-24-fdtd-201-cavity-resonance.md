# fdtd-201 cavity-resonance gate — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-24-fdtd-201-cavity-resonance-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-fdtd/tests/cavity_resonance.rs` (new) +
`crates/yee-fdtd/validation/README.md` (the `fdtd-201` row). NOTHING else.
**Out of lane** (findings, not fixes): `crates/yee-fdtd/src/**` (consume
the public API read-only — if the API can't extract a resonance, STOP +
surface), `tests/subgrid_energy_balance.rs` / the Q6 path, `fdtd-007`,
any `Cargo.toml`.

## Step ladder

### S1 — read the pattern files + API
Read `tests/fdtd_propagation.rs` (closed-PEC cavity build + field probing)
and `tests/lumped_resistor.rs` (`#[ignore]`-gated release physics-test
style). Read `src/lib.rs` for `WalkingSkeletonSolver::new` /
`step_with_source` / the field accessor used to probe E; `src/grid.rs`
`YeeGrid::vacuum`; `src/boundary.rs` `apply_pec`; `src/ntff.rs:253` for the
single-bin DFT accumulator to imitate.

### S2 — build the cavity + run + probe
New `tests/cavity_resonance.rs`: `YeeGrid::vacuum(nx,ny,nz,dx)` for a
cavity `a×b×d` with TE₁₀₁ cleanly dominant + well-separated (pick `a=d>b`;
choose `dx` so ≳15–20 cells per wavelength at f₁₀₁ for acceptable grid
dispersion). `WalkingSkeletonSolver::new`; closed PEC via `apply_pec`.
Inject an off-centre Gaussian `E_z` pulse (`step_with_source`) with σ/t0
broadband enough to cover f₁₀₁. Step N times (enough cycles for frequency
resolution — Δf ≈ 1/(N·dt)); record an interior E-field-component probe
series at a non-nodal point away from the source.

### S3 — extract the dominant resonance (no new dep)
Scan |single-bin DFT| (reuse the `ntff.rs` accumulator pattern) or an
inline Goertzel over a candidate band around f₁₀₁; peak-find the dominant
resonance. Print `extracted f / analytic f₁₀₁ / rel error`.

### S4 — assert + README
Assert the extracted dominant resonance matches analytic TE₁₀₁
`f₁₀₁=(c/2)·√((1/a)²+(1/d)²)` within the documented loose tol (≈±2–3%).
`#[ignore]`-gate the test (wall-time) with a docstring stating the tol +
the strict-±0.5%-on-refined-mesh path. Flip the `fdtd-201` README row to
live with the achieved tol.

## Verification (run in worktree; all exit 0)
```
cargo fmt --check --all
cargo clippy -p yee-fdtd --all-targets -- -D warnings
cargo test -p yee-fdtd --test cavity_resonance --release -- --ignored --nocapture   # prints + passes
cargo test -p yee-fdtd                         # rest of FDTD suite unchanged (fast, ignores the new slow test)
git diff --stat -- crates/yee-fdtd/src '**/Cargo.toml'   # must be EMPTY
```

## Escape-hatch
- If the cavity dominant-mode extraction is noisy / the wrong mode peaks,
  adjust source/probe placement + run length WITHIN this test — but if the
  public API genuinely can't probe a clean time series or the loose tol
  (±2–3%) can't be met even on a reasonable mesh, STOP and surface the
  finding (do NOT add solver `src/` code, do NOT chase ±0.5%).
- Do NOT touch the Q6 energy-balance or fdtd-007. Blocked >15 min on
  anything else → surface + stop.
- Run synchronously; no Monitor/ScheduleWakeup; no sub-agents.

## Out-of-scope (findings, not fixes)
* The strict ±0.5% refined-mesh gate (a follow-on once the loose gate lands).
* Q-factor extraction (the README row names Q-factor; the bounded first
  slice is the resonant-frequency match — note Q as a documented follow-on
  if time-boxed, do not grind on damped-mode fitting).
