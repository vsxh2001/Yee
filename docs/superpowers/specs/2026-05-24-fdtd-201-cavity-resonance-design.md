# fdtd-201 — rectangular-cavity resonance validation gate

**Status:** Draft
**Owner:** TBD
**Phase:** 2.fdtd (validation milestone `fdtd-201`)
**Type:** validation-gate addition (subsystem rotation off the MoM-port work)

## 1. Goal

Ship the `fdtd-201` validation case — the first listed Phase-2 FDTD
milestone (`crates/yee-fdtd/validation/README.md:13`: "Rectangular cavity
TE/TM Q-factor, Analytical, ±0.5%"), currently **planned but un-shipped**
(no cavity-resonance / eigenfrequency test exists in `yee-fdtd/tests/`).
Extract a rectangular PEC cavity's dominant resonant frequency from a
time-domain FDTD run and match the analytic closed-form (Pozar §6.3).
This is a published-benchmark validation case per CLAUDE.md §4 and a
genuine subsystem rotation (volumetric FDTD, not the MoM cross-section
ports the recent cycles covered).

## 2. Approach (pure consumer of the existing `yee-fdtd` public API)

1. **Cavity.** `YeeGrid::vacuum(nx, ny, nz, dx)` sized so the cavity
   `a × b × d` has a cleanly-dominant, well-separated TE₁₀₁ mode (e.g.
   `a = d > b`). `WalkingSkeletonSolver::new(grid)`; the closed PEC walls
   come from the existing outer-face clamp `boundary::apply_pec`
   (exercised in `tests/fdtd_propagation.rs`).
2. **Excitation.** An off-centre Gaussian `E_z` pulse via
   `step_with_source(i, j, k, t0, sigma)` — off-centre so it overlaps the
   TE₁₀₁ field (a centre-node placement would under-excite it).
3. **Probe.** Record a single interior E-field-component time series over
   N steps at a point away from the source (and away from mode nodes).
4. **Resonance extraction — no new dependency.** Scan the magnitude of a
   **single-bin DFT** (reuse the `ntff.rs:253` single-frequency DFT
   accumulator pattern) — or an inline Goertzel — over a candidate
   frequency band and peak-find the dominant resonance. The workspace has
   no `rustfft` (only `cufft` behind the CUDA feature), so do **not** add
   one; the single-bin-DFT scan is the sanctioned tool.
5. **Reference.** Analytic `f_{mnp} = (c/2)·√((m/a)² + (n/b)² + (p/d)²)`
   (Pozar §6.3); dominant TE₁₀₁ `f₁₀₁ = (c/2)·√((1/a)² + (1/d)²)`. Same
   closed form the FEM side already validated at fem-eig-001.

## 3. Bounded framing (tolerance policy)

**Do not promise the strict ±0.5% in the first slice.** That bound is
grid-dispersion-limited on a coarse FDTD mesh (the FEM side reached 0.5%
only after several sub-phases). Land the gate at a **documented loose
tolerance** (≈±2–3% on the dominant mode), with the strict-±0.5%-on-a-
refined-mesh path recorded in the test docstring (and, if added, an
`#[ignore]`'d strict variant). The whole test is `#[ignore]`-gated for
wall-time, mirroring the sibling slow integration tests
(`tests/lumped_resistor.rs` house style).

## 4. Definition of done

DoD-1. `crates/yee-fdtd/tests/cavity_resonance.rs`: builds the PEC
cavity, runs the pulse, extracts the dominant resonance, asserts it
matches analytic TE₁₀₁ within the documented loose tolerance. `#[ignore]`-
gated; runs green with `--release -- --ignored`.
DoD-2. The test prints a diagnostic line (extracted f vs analytic f, rel
error) so the result is visible.
DoD-3. `crates/yee-fdtd/validation/README.md` `fdtd-201` row flipped from
planned to live, citing the achieved tolerance + the refinement path.
DoD-4. No `src/` change (pure consumer of the public API — cannot regress
the solver); no new dependency; lint clean; the rest of the FDTD suite
unchanged.

## 5. NON-NEGOTIABLE

- **Tests + README only.** No `crates/yee-fdtd/src/**` edits. If the
  public API turns out to be insufficient to extract a resonance, STOP
  and surface it as a finding (do not add solver code in this track).
- Do **not** touch the deferred FDTD quagmires — the Phase 2.fdtd.7 Q6
  energy-balance (`tests/subgrid_energy_balance.rs`, 75–79% drift) and
  `fdtd-007` (wrong-reference, `validation/README.md:42`). This track is
  independent of both.
- No new `Cargo.toml` dependency (no `rustfft`); reuse the single-bin DFT.

## 6. References

* Pozar, *Microwave Engineering* §6.3 (rectangular cavity resonances).
* `crates/yee-fdtd/validation/README.md:13` (the `fdtd-201` contract).
* Pattern files: `crates/yee-fdtd/tests/fdtd_propagation.rs` (closed-PEC
  cavity setup), `crates/yee-fdtd/tests/lumped_resistor.rs`
  (`#[ignore]`-gated release physics-test house style).
* API: `crates/yee-fdtd/src/lib.rs` (`WalkingSkeletonSolver::new`,
  `step_with_source`), `src/grid.rs` (`YeeGrid::vacuum`), `src/boundary.rs`
  (`apply_pec`), `src/ntff.rs:253` (single-bin DFT to reuse).
