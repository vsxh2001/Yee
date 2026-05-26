# ADR-0073 — Phase 2.fdtd.py.1: FDTD Cavity-Resonance Python Driver

**Date:** 2026-05-26
**Status:** Accepted
**Phase:** 2.fdtd.py.1

---

## Context

The fdtd-201 (TE₁₀₁ rectangular cavity resonance, ADR-0062) and fdtd-201.x (TE₂₀₁
higher-mode, ADR-0066) gates are `#[ignore]`-gated Rust integration tests that validate the
FDTD solver's ability to extract cavity resonant frequencies via a broadband DFT scan.
Neither gate is:

- callable from Python (unlike fdtd-202 / `run_cavity_q` shipped in ADR-0072), or
- registered in `Report::run_all()` (so `yee validate all` does not surface them).

The pattern established by Phase 2.fdtd.py.0 (ADR-0072) is: expose a self-contained Python
callable in `yee-py/src/fdtd.rs` that re-expresses the validation gate logic, then register
the Rust gate as Skipped in `run_all`.

---

## Decision

Implement **Phase 2.fdtd.py.1** with the following artifacts:

1. **`run_cavity_resonance()` / `PyCavityResonanceResult`** in `crates/yee-py/src/fdtd.rs` —
   a Python-callable FDTD rectangular cavity simulation that extracts TE₁₀₁ resonant
   frequency via DFT scan and verifies within ±2.5 % of the analytic Pozar §6.3 value.

2. **Aggregator registrations** in `crates/yee-validation/src/lib.rs`:
   - `run_fdtd_201_cavity_resonance()` → `Skipped` (wall-time ~5–15 s release)
   - `run_fdtd_201x_cavity_higher_mode()` → `Skipped` (same wall-time class)

3. **Tutorial 11** at `docs/src/tutorials/11-fdtd-cavity-resonance-from-python.md`.

No `yee-fdtd/src/` changes; the existing `#[ignore]`-gated Rust tests are untouched.

---

## Alternatives considered

**A. Move logic to yee-fdtd/src/ and expose through pub API.** Cleaner in principle, but
adds a `run_cavity_resonance` public function to a crate whose scope is the solver building
blocks — not simulation drivers. The fdtd.py.0 pattern (logic in yee-py) keeps driver code
in the bindings layer where it belongs. Rejected.

**B. Un-ignore the Rust tests (register as non-Skipped in run_all).** The 5–15 s wall-time
makes the default `yee validate all` ~2× slower. fdtd-202 at 0.38 s was the threshold; we
keep both fdtd-201 and fdtd-201.x as Skipped until a faster reduced-grid variant justifies
promoting them. Rejected for now.

**C. Add fdtd-201.x Python driver in this increment.** Scope creep; the TE₂₀₁ higher-mode
test adds little additional physics coverage at this stage. Deferred to Phase 2.fdtd.py.2.

---

## Consequences

- `from yee import run_cavity_resonance, CavityResonanceResult` works after `maturin develop`.
- `yee validate all` output gains two `SKIP` lines for fdtd-201 and fdtd-201.x, making the
  aggregator's FDTD coverage complete through the shipped gates.
- Tutorial 11 gives Python users an interactive way to explore FDTD cavity physics.
- No existing gates are weakened; no existing tests change.
