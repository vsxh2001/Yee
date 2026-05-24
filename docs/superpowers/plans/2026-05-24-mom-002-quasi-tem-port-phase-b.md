# mom-002 quasi-TEM port Phase B — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-24-mom-002-quasi-tem-port-phase-b-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-mom/src/ports.rs` (the `NumericalCrossSection`
quasi-TEM solve path + mode_profile caching), `crates/yee-mom/tests/mom_002_numerical_waveport.rs`
(extend to Part B), `crates/yee-mom/src/lib.rs` (only a test-only
`__internal` helper if needed), `ROADMAP.md`, `docs/src/decisions/0061-*.md`.
**Out of lane** (findings, not fixes): the eigensolver
(`solve.rs` `solve_dense_mixed`/`_quasi_tem`/`cutoff_candidates`,
`reference.rs`, the element matrices — consume read-only), the mom-002
kernel/Greens/gate/tripwire, `crates/yee-validation/src/lib.rs` solve
internals (constants read-only), `crates/yee-fem/**`, `crates/yee-py/**`.
No `Cargo.toml` dependency.

## Step ladder

### A — NumericalCrossSection quasi-TEM solve path (production wiring)
`ports.rs`: add a quasi-TEM solve path — mirror `solve`'s First branch
(scatter the eigenvector → `mode_profile`, build `tri_edges_cache`) but
call `solve_dense_mixed_quasi_tem` (which returns the same
`MixedEigenSolution`). Expose via `solve_quasi_tem(freq)` or a mode
selector; First-order closed-guide `solve` stays default + bit-identical.
Smoke: the microstrip cross-section solves quasi-TEM + `mode_profile`
populated + `e_tangential_at` returns finite E_t.

**Verification:** `cargo test -p yee-mom --lib ports` (or the smoke) green;
closed-guide `solve` unchanged.

### B — mom-002 numerical port + Z_in (HARD ~30-min coupling cap)
Extend `mom_002_numerical_waveport.rs`: build the mom-002 microstrip
cross-section, `solve_quasi_tem`, feed `mode_profile` via
`WavePort::with_numerical_cross_section` + Numerical2D RHS to the mom-002
line, LU-solve, extract `|Z_in|`. **Cap:** if the cross-section→RWG modal
RHS for a microstrip needs glue that doesn't exist → document the blocker
+ STOP. Report `|Z_in|` vs 674 Ω vs 51 Ω as a NON-FAILING diagnostic
(no mom-002 re-gate).

**Verification:** the diagnostic prints the comparison (or the blocker
finding); mom-002 gate + mom-001/003 unchanged.

### C — ADR-0061 + ROADMAP
Record the outcome + recommendation (adopt / residual-not-port /
glue-needed).

## Full verification (all exit 0)
```
cargo fmt --check --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p yee-mom --test mom_002_numerical_waveport
cargo test -p yee-mom --test eigensolver_microstrip_quasi_tem   # quasi-TEM gate unchanged
cargo test -p yee-mom --test eigensolver_wr90 --test te10_waveport --test wave_port_numerical_te10  # closed-guide + port path
cargo test -p yee-mom --lib eigensolver                         # eigensolver unchanged
git diff --stat -- '**/Cargo.toml'        # expect EMPTY
```
(mom-001/002 full gates slow — rely on CI for the full mom suite.)

## Escape-hatch
Part B's hard ~30-min coupling cap is the escape-hatch: wired → the Z_in
comparison; not wired → documented glue-blocker + stop. If Part B's Z_in
points back at the kernel/Greens (not the port), STOP + document — do
NOT re-open the forensics. If the numerical port's Z_in is no better than
the delta-gap, that is a valid finding (the residual is not the port) —
do NOT chase. NEVER touch the eigensolver internals / the mom-002 gate /
the verified reference.

## Out-of-scope (findings, not fixes)
* Adopting the numerical port as the mom-002 PRODUCTION excitation +
  re-gating mom-002 (a separate track if Part B shows a clear win).
* The mom-002 kernel/Greens. yee-fem; yee-py.
