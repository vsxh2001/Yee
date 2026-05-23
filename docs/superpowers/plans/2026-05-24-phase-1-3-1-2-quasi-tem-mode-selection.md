# Phase 1.3.1.2 — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-24-phase-1-3-1-2-quasi-tem-mode-selection-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-mom/src/eigensolver/{solve,mod}.rs`,
`crates/yee-mom/tests/eigensolver_microstrip_quasi_tem.rs` (new),
`crates/yee-mom/src/eigensolver/reference.rs` (ONLY to add a closed-form
HJ ε_eff helper if needed — do NOT touch the verified slab-loaded
transcendental), `ROADMAP.md`, `docs/src/decisions/0060-*.md`.
**Out of lane** (findings, not fixes): the slab-loaded reference
dispersion (verified — do not alter), `assemble_mixed`/`assemble_mixed_p2`
element matrices (consume), `ports.rs` public contract, `crates/yee-fem/**`,
`crates/yee-py/**`. No `Cargo.toml` dependency.

## Step ladder (FEASIBILITY-FIRST)

### Q1 — FEASIBILITY (bounded ~40 min): surface the quasi-TEM mode
Build a (shielded) canonical microstrip cross-section (strip width w,
substrate height h, ε_r, large air box, PEC outer). Extend
`cutoff_candidates` / the selection so a `k_c²≈0` transverse-dominated
PROPAGATING candidate is gathered (near-zero shift-invert rung and/or a
TEM-like seed) and survives the converged-eigenvector transverse screen
(the gradient nulls, E_t≈0, are screened out). **DECISION (bounded):**
if a transverse-dominated quasi-TEM mode (k_c²≈0, β>0, ε_eff in the
microstrip ballpark) is surfaced → Q2. If not within the budget →
document the specific gathering blocker, STOP, queue a follow-on. Do NOT
grind into a multi-step chase.

**Verification:** the microstrip solve returns a transverse-dominated
β>0 quasi-TEM mode (not "no propagating cutoff candidate").

### Q2 — HJ validation (DoD-2)
`eigensolver_microstrip_quasi_tem.rs`: assert the quasi-TEM
`ε_eff=(β/k₀)²` matches the Hammerstad-Jensen `ε_eff` (canonical
microstrip closed form) within ≤5–10% (loose, box-truncation-perturbed).
Document the box size + tolerance rationale.

**Verification:** the microstrip-vs-HJ test passes within tol.

### Q3 — No-regression guard (DoD-3)
WR-90 TE10, FR-4 §4 gate, homogeneous canary, uniform anchor,
vertical-slab, coupling guards — bit-identical / unchanged. If a unified
selection relaxation regresses any, scope the quasi-TEM path to a
separate entry-point keeping the closed-guide selection bit-identical.

**Verification:** the full block below, all exit 0.

### Q4 — ROADMAP + ADR-0060

## Full verification (all exit 0)
```
cargo fmt --check --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p yee-mom --lib eigensolver
cargo test -p yee-mom --test eigensolver_microstrip_quasi_tem   # new
cargo test -p yee-mom --test eigensolver_inhomogeneous
cargo test -p yee-mom --test eigensolver_wr90
cargo test -p yee-mom --test wave_port_numerical_te10 --test te10_waveport
git diff --stat -- '**/Cargo.toml'        # expect EMPTY
```

## Escape-hatch
The Q1 bounded feasibility cap IS the escape-hatch: surfaced → Q2; not
surfaced → documented gathering-blocker finding + stop (a follow-on
track scopes the harder near-zero separation). If surfaced but HJ is
off by ≫10% (not box-truncation-explicable), document — do NOT chase it
into a multi-step refinement (that would be the ε_r=10.2-style trap).
NEVER weaken WR-90 / FR-4 / homogeneous; NEVER alter the verified
slab-loaded reference.

## Out-of-scope (findings, not fixes)
* Adopting the quasi-TEM port into the mom-002 production excitation (a
  follow-on once the quasi-TEM mode is validated).
* True open-boundary (PML/absorbing) microstrip (vs the shielded box).
* yee-fem; yee-py.
