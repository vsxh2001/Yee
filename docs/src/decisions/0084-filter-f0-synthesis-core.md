# ADR-0084: Filter Phase F0 — synthesis core (`yee-synth` + `yee-filter`)

**Status:** Accepted
**Date:** 2026-05-29
**Supersedes:** none
**Related:** `FILTER-DESIGN-ROADMAP.md` (Phase F0); ADR-0028/0031 (NL design
surface, `yee-design` — distinct: NL→intent, not filter synthesis)

---

## Context

`FILTER-DESIGN-ROADMAP.md` sets the project's final goal: end-to-end RF filter
design. Its first phase, **F0**, is the synthesis *walking skeleton* — the
minimal end-to-end pipe `FilterSpec → ideal response`, pure math, no EM, no
layout, no new heavy dependency. Everything downstream (planar/waveguide/lumped
back-ends, surrogate dimensional synthesis, layout, export) plugs into the data
model and the synthesis output that F0 establishes.

There is no filter-synthesis code in the tree today; `yee-design` is the NL
intent surface (emit/estimate/intent/offline), unrelated.

## Decision

Create two new pure-math crates and wire a CLI entry:

- **`yee-synth`** — classical filter synthesis, no EM, no I/O:
  - Lowpass-prototype g-values: **Butterworth** (`g_k = 2·sin((2k−1)π/2N)`,
    `g0=g_{N+1}=1`) and **Chebyshev** (Pozar §8.3 eq 8.53 recursion; see spec).
  - Lowpass→bandpass frequency transform.
  - All-pole **coupling matrix** + external-Q synthesis from g-values
    (`k_{i,i+1}=FBW/√(g_i g_{i+1})`, `Qe1=g0 g1/FBW`, `Qen=g_N g_{N+1}/FBW`).
- **`yee-filter`** — filter-domain data model + ideal response + flow scaffold:
  - Types: `FilterSpec`, `Approximation`, `Prototype`, `CouplingMatrix`,
    `Topology`, `FilterProject` (all `serde`).
  - Ideal bandpass response from the synthesized prototype (closed-form
    Chebyshev/Butterworth transfer function) → S-parameter sweep; spec-mask
    pass/fail (`SpecMask`: passband ripple/RL, stopband rejection points).
- **`yee-cli`**: `yee filter synth <spec.toml>` → prints the prototype +
  coupling matrix, writes the ideal S-params as Touchstone (via `yee-io`), and
  reports spec-mask pass/fail. Exit 0 on success.

Workspace `Cargo.toml` gains the two members. Only `nalgebra` (already in tree)
is needed for linear algebra; **no new dependency**.

### Why closed-form response (not coupling-matrix→S) for F0's `filt-001`

The coupling-matrix → S-parameter path (Hong-Lancaster `[A]=[q]+pU−jM`) is the
general route and lands in a later phase. For the F0 skeleton, the synthesized
prototype's **closed-form** transfer function (`|S21|²=1/(1+ε²T_N²(Ω))`) is
exact, trivially validated, and sufficient to prove the `spec → response → mask`
pipe. The coupling matrix is still synthesized and emitted (gate `synth-002`);
driving S-params *from* it is F1+ work.

## Consequences

**Ships:** `yee-synth`, `yee-filter`, `yee filter synth` CLI, and three
published-benchmark gates (see spec §DoD): `synth-001` (Butterworth + Chebyshev
g-values vs Pozar Tables 8.3/8.4 / Matthaei-Young-Jones, ≤1e-3), `synth-002`
(all-pole coupling coefficients + Qe vs a worked Hong-Lancaster example),
`filt-001` (synthesized Chebyshev ideal response meets its own ripple/RL/
rejection mask). These register in `yee-validation` as fast (`Run`) cases.

**Not in scope (later phases):** coupling-matrix→S-parameter realization,
elliptic/cross-coupled (Cameron) synthesis, any EM, any layout/export, the GUI
wizard. F0 is math + data model + CLI only.

**No new dependency.** Lane: `crates/yee-synth/**`, `crates/yee-filter/**`,
workspace `Cargo.toml`, `crates/yee-cli/**`, and a `synth-*`/`filt-*`
registration in `crates/yee-validation/**`.

---

## References

- Pozar, *Microwave Engineering* 4e, §8.3–8.4 (prototype g-values, Tables 8.3/8.4).
- Matthaei, Young & Jones, *Microwave Filters, Impedance-Matching Networks…*,
  Table 4.05-2 (Chebyshev g-values).
- Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications*, ch. 8
  (coupling coefficients, external Q).
- `FILTER-DESIGN-ROADMAP.md` (Phase F0); `docs/superpowers/specs/2026-05-29-filter-f0-synthesis-core-design.md`;
  `docs/superpowers/plans/2026-05-29-filter-f0-synthesis-core.md`.
