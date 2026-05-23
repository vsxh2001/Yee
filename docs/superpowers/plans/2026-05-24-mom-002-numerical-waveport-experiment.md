# mom-002 numerical-microstrip-wave-port — experiment plan

**Spec:** `docs/superpowers/specs/2026-05-24-mom-002-numerical-waveport-experiment-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-mom/src/ports.rs` (ONLY if the `Numerical2D`
microstrip-into-planar-MoM coupling needs a small, documented helper —
prefer a test-side construction first), `crates/yee-mom/tests/` (a new
`mom_002_numerical_waveport.rs` diagnostic), `crates/yee-validation/`
(only to read the mom-002 constants / `z_in_with_greens_tem` — read-only
unless a non-failing diagnostic is added there), `docs/src/decisions/0059-*.md`,
`ROADMAP.md`.
**Out of lane** (findings, not fixes): the cross-section eigensolver
(`reference.rs`, `assembly.rs`, `solve.rs` — read-only), the mom-002
kernel / Greens / gate / tripwire band (do NOT touch), `crates/yee-fem/**`,
`crates/yee-py/**`. No `Cargo.toml` dependency.

## Step ladder

### A — Feasibility (HARD 30-min cap)
1. Build a microstrip cross-section `TriMesh2D` (FR-4 substrate ε_r≈4.3 +
   signal strip + ground + air box, the transverse plane of the mom-002
   line; strip width 2.94 mm per `MOM_002_STRIP_WIDTH_M`). Tag materials.
2. `NumericalCrossSection::new(...).solve(freq)` — confirm a physical
   quasi-TEM β / Z_w (sanity: ε_eff in the FR-4 ballpark ≈3.3).
3. Attempt to feed its modal `E_t` to the mom-002 MoM line via
   `WavePort::with_numerical_cross_section` + `ModalDistribution::Numerical2D`
   → a port RHS for the planar-MoM solve.
4. **DECISION (hard cap):** if the 2-D-cross-section→2.5-D-RWG-port
   coupling can be wired cleanly → Phase B. If not within 30 min →
   STOP, write the finding (what glue the `Numerical2D` arm lacks for
   microstrip ports), commit it, report. Do NOT force the coupling.

**Verification (A):** the cross-section solve returns a physical mode;
the feasibility decision is recorded.

### B — Comparison (only if A wired)
1. Solve the mom-002 line with the numerical wave-port RHS; extract
   `|Z_in|`. Compare to the delta-gap baseline 674 Ω + HJ target ≈51 Ω.
2. Ship as a NON-FAILING diagnostic test (`mom_002_numerical_waveport.rs`)
   — print numerical-port `|Z_in|`, 674, 51. Do NOT assert a tripwire /
   do NOT touch the mom-002 gate.

**Verification (B):** the diagnostic runs + prints the comparison;
mom-002 gate + mom-001/003 unchanged.

### C — ADR-0059 + ROADMAP
Record the outcome + recommendation (adopt-numerical-port follow-on /
residual-not-the-port / port-infra-glue-needed).

## Full verification (all exit 0)
```
cargo fmt --check --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p yee-mom --test mom_002_numerical_waveport   # the new diagnostic (or feasibility finding)
cargo test -p yee-mom --test te10_waveport --test wave_port_numerical_te10  # no regression to the port path
git diff --stat -- '**/Cargo.toml'        # expect EMPTY
```
(mom-001/002/003 full gates are slow — run the targeted port + the new
diagnostic; rely on CI for the full mom suite.)

## Escape-hatch
The HARD 30-min Phase-A cap IS the escape-hatch: a clean wiring → B; no
clean wiring → documented finding + stop. Beyond that: if Phase B's
numerical-port `|Z_in|` points back at the kernel/Greens (not the port),
STOP + document — do NOT re-open the forensic analysis. NEVER touch the
eigensolver, the mom-002 gate, or the kernel.

## Out-of-scope (findings, not fixes)
* Adopting the numerical port as the mom-002 production excitation (a
  follow-on track if the experiment shows a clear win).
* The mom-002 kernel / Greens / tripwire band. yee-fem; yee-py.
