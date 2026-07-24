# ADR-0226: FS.4.2b — per-layer stackup loss: σ map + stripline attenuation gate

**Date:** 2026-07-23 · **Status:** accepted · **Track:** FS.4 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-23-fs4-2b-stackup-loss-design.md`
**Plan:** `docs/superpowers/plans/2026-07-23-fs4-2b-stackup-loss.md`
**Predecessors:** FS.4.0 (ADR-0215) — `Stackup` + `voxelize_stackup`; FS.4.1 (ADR-0221) —
vias; FS.4.2a (ADR-0225) — H-probes + `engine-stripline-z0-001`; FS.2c —
`substrate_sigma_cells` (single-layer tan δ over the whole substrate ε-map).

## Context

`StackupLayer.loss_tangent` has existed since FS.4.0 but nothing consumed it —
`substrate_sigma_cells` (FS.2c) maps exactly one tan δ over an entire single-layer
substrate. Multilayer boards mix materials (an FR-4 core at tan δ ≈ 0.02 under a
low-loss prepreg at tan δ ≈ 0.004 is a routine real stackup); per-layer dielectric
loss is table stakes for FS.4's "every real board" story, and every downstream FS.4
gate to date (`engine-stripline-eeff-001`, `engine-stripline-z0-001`,
`engine-stackup-via-001`) has been lossless.

## Decision

### 1. `stackup_sigma_cells` (Task 1, commit `54a974e`)

`pub fn stackup_sigma_cells(model: &MicrostripModel, stackup: &Stackup, f_ref_hz: f64)
-> Vec<f64>` in `yee-voxel`: per-cell σ = `2π f_ref ε₀ ε_r(layer) tan δ(layer)`,
assigned by re-deriving each cell's k-band from `stackup`'s layer heights and the
model's `dz` — the same bookkeeping `voxelize_stackup` uses to fill ε — rather than
inferring the layer from the already-written ε value (ε values can coincide across
layers, so that inference would be ambiguous). Air and metal (trace/ground/lid)
cells get σ = 0.

**A real bug caught by the consistency test, not assumed away**: an uncommitted
first-draft implementation (found already sitting in the working tree at the start
of Task 1, no test coverage) assigned σ by `k` alone, ignoring that the ε array's
`nx+1`/`ny+1` shape carries one **padding plane** per axis that `voxelize_stackup`'s
and `voxelize_microstrip`'s fill loops never write (only `i in 0..nx, j in 0..ny`) —
it stays at the ε = 1.0 (air) default. The draft leaked the layer's σ onto that
padding plane at every k inside a filled band, a real footprint bug past where ε was
ever actually set. `single_layer_stackup_sigma_matches_substrate_sigma_cells` caught
it directly: `substrate_sigma_cells` derives σ from the *actual* ε value (correctly
0 at the padding plane), so the two vectors diverged at the first padding-plane cell.
Fixed by mirroring the exact `ii < nx && jj < ny` footprint check the fill loops use
before the k-band lookup, else 0.0 — the docstring's claim of "0 elsewhere (air, the
lid plane, any plane above the stack)" is now actually true end-to-end.

**Unit tests** (`rf_tool_tests`, reusing the `voxel_stackup_002.rs` two-layer fixture
idiom):

- `stackup_sigma_cells_matches_each_layer_band_exactly` — two-layer stack (ε_r
  2.2/tan δ 0.02, 2 cells; ε_r 4.4/tan δ 0.005, 3 cells) → exact σ per band,
  boundary k exact (not blended), 0 above the stack.
- `all_zero_loss_tangent_is_a_provable_no_op` — same geometry, tan δ = 0 both
  layers → every returned σ is exactly `0.0`. This is the loss-off no-op guarantee
  the spec requires, proven at the unit level (cheap) rather than only inferred from
  the full-wave gate.
- `single_layer_stackup_sigma_matches_substrate_sigma_cells` — single-layer
  `Stackup` vs the equivalent `Substrate` → **bit-identical** vectors (both compute
  the same `omega * EPS0 * eps_r * tan_d` expression, and a single-layer
  `voxelize_stackup` model is bit-identical to `voxelize_microstrip` per
  `voxel-stackup-001`).

No conductor loss is modeled (PEC everywhere, per the spec's non-goals).

### 2. Gate `engine-stripline-alpha-001` (Task 2, commit `f3f1ed3`)

New file `crates/yee-engine/tests/stripline_alpha.rs`. Reuses
`engine-stripline-z0-001`'s exact cross-section (ε_r 2.2, b = 16 cells = 3.2 mm, w/b
= 0.8125, dx = 0.2 mm) — `Stackup::symmetric_stripline` hardcodes `loss_tangent:
0.0`, so the gate hand-builds a two-half-layer `Stackup` with tan δ = 0.02 and feeds
it through `stackup_sigma_cells(&model, &stack, F_REF_HZ)` into
`Materials::sigma_cells`, the same E.1 lossy CA/CB update path `engine-loss-001`
already exercises for microstrip.

**Closed form**: stripline is pure TEM entirely inside the dielectric (no fringing
into air, unlike microstrip), so the attenuation closed form is exact — not an
approximation like Pozar §3.199's ε_eff-fudged microstrip form:

```text
α_d = (π f √ε_r / c) · tan δ   [Np/m]      (× 8.686 → dB/m)
```

**Extraction**: two `Ez`-column measurement planes (the `engine-stripline-z0-001`
V-column idiom — ground k=0 up to, excluding, the trace plane, summed × dz) on one
FDTD run, plane A at 2.0 λg and plane B at 4.0 λg downstream of the port.
`α_meas = ln(|V_A|/|V_B|) / (x_B − x_A)`, `x_B − x_A` from the actual
grid-quantized index separation. Both planes see the identical launched wave on a
single pass — this is launch-normalized by construction (one run, one wave, two
taps), so ADR-0204's warning against absolute single ratios *across separate runs*
with different incident waves does not apply here; there is no second run to
diverge from. Gated at 0.9× the wall-reflection time computed from the farther
plane, before the hard-PEC end-wall reflection reaches either plane. Line length
bumped to 8.5 λg (vs Z₀'s 8 λg) to leave enough margin for plane B's Gaussian-pulse
tail to clear before the gate closes — sized analytically before running any FDTD,
per the `stripline_eeff.rs` "window hygiene" lesson (a too-short margin there
inflated a measurement by 14.5%).

**Lossless control, same fixture**: tan δ = 0.0 through the identical
`stackup_sigma_cells` path, in the same test — the differential that kills
systematic gating bias. Asserts the measured α stays under 5% of the lossy run's
α_ref (the no-op bound) *and* that plane A's lossless `|V|` phasor lands in the
same sane, non-trivial range `engine-stripline-z0-001` measures (~2.33) — proving
σ = 0 didn't silently zero the field rather than genuinely not attenuate it.

**Constant-σ vs true tan δ (documented, not asserted on)**: `stackup_sigma_cells`
maps tan δ to σ at a single reference frequency
(`σ = 2π f_ref ε₀ ε_r tan δ`); the FDTD update then treats σ as
frequency-*independent*, so the discrete model's implied loss tangent drifts as
`tan δ_eff(f) = tan δ(f_ref) · f_ref / f` away from f_ref — a real modeling
deviation off-reference, not a bug (true Debye/dispersive tan δ is a separate,
unshipped lane — the dispersive ADE materials track already covers frequency-
dependent loss where that fidelity is needed). The gate grades only at f_ref, where
the constant-σ model is exact; it does not sweep frequency or assert on the
deviation.

## Measured result

```
$ cargo test -p yee-engine --release --test stripline_alpha -- --ignored --nocapture
engine-stripline-alpha-001: tan_d = 0.02, d = 67.40 mm | |V_A| = 2.0486e0, |V_B| = 1.8002e0
  -> alpha_meas = 1.9178 Np/m (16.6581 dB/m) vs closed form alpha_ref = 1.8652 Np/m
     (16.2010 dB/m) -> err 2.821 %
  lossless control: |V_A| = 2.3343e0, |V_B| = 2.3334e0 -> alpha_lossless = 0.005630 Np/m
     (0.0489 dB/m) — numeric gating-bias floor
test stripline_alpha_matches_the_pozar_dielectric_loss_closed_form ... ok   (115.68 s)
```

**2.821% measured**, well inside the plan's ≤ 10% target and the ≤ 20%
STOP-and-root-cause threshold — no widening needed. The lossless floor
(0.00563 Np/m, 0.30% of α_ref) is well inside the 5% no-op bound; plane-A's
lossless `|V|` = 2.33 matches `engine-stripline-z0-001`'s measured `|V| ≈ 2.33`
almost exactly (same fixture geometry and drive amplitude, as expected).

## Tolerance pinned

`rel_err ≤ 0.10` (10%, matches the plan's explicit target) on α at f_ref; lossless
no-op bound ≤ 5% of α_ref. Measured 2.821% / 0.30% give substantial headroom on
both. Never widened — not needed, since the first measurement already cleared
tolerance.

## Bit-exactness / regression discipline (unmodified gates, every commit)

The binding gate command —
`cargo test -p yee-compute --release --test graded_uniform_bitexact --test
gpu_graded_parity --test gpu_cpu_parity -- --include-ignored` — stayed green after
both Task 1 and Task 2 commits (real GPU adapter `NVIDIA GeForce RTX 5060 Ti`, not
SKIPPED). `cargo test -p yee-voxel --release` green throughout, including Task 1's
3 new unit tests. Both stripline full-wave gates re-run unmodified after Task 2:
`engine-stripline-z0-001` (1.271%, unchanged) and `engine-stripline-eeff-001`
(0.065%, unchanged) — `git diff` on `stripline_z0.rs`/`stripline_eeff.rs` empty
throughout. Workspace clippy (default + `--no-default-features` on `yee-compute`)
and `cargo fmt --check --all` clean before every commit; `missing_docs` clean.

## Runtime note

The gate measured 115.68 s (two full FDTD runs — lossy + lossless — in one test
function, vs the Z₀ gate's single run at ~24-40 s), above the plan's "~2×" runtime
guideline. Root cause (analyzed, not hand-waved): the differential design
inherently needs two runs, and the pulse-decay margin required to avoid the
`stripline_eeff.rs`-documented truncation-bias failure mode forced the longer 8.5 λg
box. Smaller geometries were evaluated analytically (down to ~20% cheaper) but every
one cut the margin below the validated ~11% safety margin used here; given "never
weaken any assertion" / "root-cause, never widen", correctness margin was kept over
shaving CI seconds. This is a small fraction of the yee-engine blanket
`--include-ignored` job's total budget and required no `ci.yml` change — the gate's
name doesn't match any of the job's `--skip` prefixes (`antenna_`, `patch_`,
`inset_`, `design_loop_`, `graded_`), confirmed by grep, so it is auto-picked-up the
same way `stripline_z0`/`stripline_eeff` were.

## Verdict

**GO.** `stackup_sigma_cells` closes the FS.4 per-layer-loss gap with a provable
loss-off no-op and a bit-identical single-layer consistency pin against the
pre-existing FS.2c path; `engine-stripline-alpha-001` confirms the σ-map is
physically correct end-to-end on the first honest full-wave measurement (2.821%
against an exact TEM closed form, no root-cause detour needed), on the exact same
fixture the Z₀ and ε_eff gates already validate. FS.4's stripline trio
(ε_eff / Z₀ / α) is now complete against three independent exact closed forms.

## What remains of FS.4.2 (queued, not attempted here)

Per the spec's non-goals: frequency-dependent σ(f) / Debye fitting (the dispersive
ADE materials lane already exists for that fidelity, not re-derived here); conductor
(finite-σ metal) loss; a MoM cross-check of the dielectric-loss model; automesh
awareness of per-layer loss tangent in the stackup rulebook. None of these were
scoped into FS.4.2b.
