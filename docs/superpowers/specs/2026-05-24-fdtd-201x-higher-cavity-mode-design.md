# fdtd-201.x — higher-order rectangular-cavity mode resonance gate

**Status:** Draft
**Owner:** TBD
**Phase:** 2.fdtd validation (sibling to fdtd-201, ADR-0062)
**Type:** validation-gate addition (tests-only; near-zero-risk reuse of the
fdtd-201 harness)

## 1. Goal

fdtd-201 validated the **dominant** TE₁₀₁ cavity mode. Add a gate for a
**higher-order** mode (TE₂₀₁) extracted from the same time-domain FDTD
cavity harness, validating mode **selectivity** + the solver's
higher-frequency grid-dispersion behaviour — a distinct claim from "the
dominant mode is right". Reference: analytic Pozar §6.3 (already coded at
`crates/yee-fdtd/tests/cavity_resonance.rs:167`).

## 2. Approach (clone the fdtd-201 harness; tests-only)

Reuse `cavity_resonance.rs`'s machinery verbatim: `YeeGrid::vacuum` +
`WalkingSkeletonSolver::new` + `boundary::apply_pec` (closed PEC cavity),
off-centre Gaussian into `grid.ey` via the custom step body, interior
E_y probe, single-bin-DFT frequency scan, peak-find.

**Break the degeneracy with a ≠ d.** For a cubic-in-x/z box (a = d),
TE₂₀₁ and TE₁₀₂ are degenerate (both `√5·f₁₀₁`) — a peak can't be
attributed to a *named* mode. Use **a ≠ d** (e.g. a 24×10×16 grid →
a = 0.24, b = 0.10, d = 0.16 m at dx = 10 mm; cells stay cubic, only the
counts differ) so the target mode is cleanly isolated in frequency.
Pick the target (TE₂₀₁) + place the source/probe at one of its E_y
antinodes (TE₂₀₁: `E_y ∝ sin(2πx/a)·sin(πz/d)`, antinodes at x = a/4, 3a/4)
so it couples; widen the DFT scan band to bracket the target mode and
exclude the dominant TE₁₀₁ + neighbours. Verify the target is cleanly
separated (no other mode within the scan band).

## 3. Definition of done

DoD-1. New gate (`crates/yee-fdtd/tests/cavity_higher_mode.rs`, or a sibling
`#[test] #[ignore]` fn in `cavity_resonance.rs`) extracts the chosen
higher-order mode's resonance + asserts it matches the analytic Pozar §6.3
`f_mnp` within a documented loose tolerance (≈±2.5%; grid dispersion is the
floor + is worse at higher f, so do NOT promise ±0.5%). Prints the
extracted-vs-analytic diagnostic line.
DoD-2. The geometry makes the target mode non-degenerate + cleanly
separated in the scan band (a ≠ d); the gate names the specific mode it
validates.
DoD-3. `#[ignore]`-gated (wall-time); `crates/yee-fdtd/validation/README.md`
gets a row for the higher-mode gate.
DoD-4. No `src/` change (pure consumer of the existing public API); no new
dependency; lint clean; the fdtd-201 gate + rest of the suite unchanged.

## 4. NON-NEGOTIABLE

- **Tests + README only.** No `crates/yee-fdtd/src/**` edits. If the public
  API can't isolate/extract the higher mode, STOP + surface (don't add a
  mode-decomposition probe to `src/`).
- Do NOT touch the quagmires: the subgrid surface (`subgrid_*`,
  `berenger_*`, Q6 energy-balance) and fdtd-007.
- No new `Cargo.toml` dependency.

## 5. References

* Pozar §6.3 (`f_mnp = (c/2)√((m/a)²+(n/b)²+(p/d)²)`).
* Pattern: `crates/yee-fdtd/tests/cavity_resonance.rs` (fdtd-201 — the exact
  harness to clone), `tests/lumped_resistor.rs` (`#[ignore]` release-physics
  style). `crates/yee-fdtd/validation/README.md` (add the row).
* API: `YeeGrid::vacuum` (`src/grid.rs`), `WalkingSkeletonSolver` step
  primitives (`src/lib.rs`), `boundary::apply_pec` (`src/boundary.rs`).
