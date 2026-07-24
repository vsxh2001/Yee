# ADR-0227: FS.4.2c — automesh stackup integration: N-layer rulebook + lid b/16 rule

**Date:** 2026-07-24 · **Status:** accepted · **Track:** FS.4 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-24-fs4-2c-automesh-stackup-design.md`
**Plan:** `docs/superpowers/plans/2026-07-24-fs4-2c-automesh-stackup.md`
**Predecessors:** FS.0a `auto_dx` rulebook (ADR-0204); FS.4.0 — `Stackup` +
`voxelize_stackup` + the measured lesson "confined lidded modes need ≥ ~16 cells
across b, else β reads 7 % high" (ADR-0215); FS.4.2a/b — stripline Z₀/α gates
(ADR-0225/0226).

## Context

`yee_engine::automesh::auto_dx(layout, f_max)` knows one substrate
(`layout.substrate`). Stackup boards (FS.4.0) have N layers + an optional lid, and
the ADR-0215 b-resolution lesson ("lidded/confined modes need ≥ ~16 cells across
the confinement dimension") lived only in prose plus a hand-set fixture constant —
every stackup full-wave gate to date (`engine-stripline-eeff-001`,
`engine-stripline-z0-001`, `engine-stripline-alpha-001`) hand-picked `dx = 0.2 mm`
rather than deriving it. Push-button meshing (the FS.0 wedge) needs to extend to
multilayer before the rule can be trusted on a board this repo hasn't hand-tuned yet.

## Decision

### 1. `auto_dx_stackup` (Task 1, commit `c88149c`)

`pub fn auto_dx_stackup(layout: &Layout, stackup: &yee_layout::Stackup, f_max_hz: f64)
-> f64` in `yee-engine`, immediately after `auto_dx`. Largest dx satisfying, same
`[1 µm, 1 mm]` clamp as `auto_dx`:

- **wavelength**: `dx ≤ λ_min/20`, `λ_min = c/(f_max·√ε_r_max)`, `ε_r_max` = the max
  `eps_r` over `stackup.layers` (generalizes `auto_dx`'s single-substrate term).
- **per-layer resolution**: `dx ≤ h_i/3` for **every** layer, fold-min (generalizes
  the h/3 rule — a buried interface under-resolved is the same silent failure as an
  under-resolved single substrate).
- **feature**: `dx ≤ min_feature_m(layout)/2`, reusing `min_feature_m` unchanged.
- **lid rule (the ADR-0215 lesson made a rule)**: if `stackup.lid`, additionally
  `dx ≤ b/16` where `b = stackup.total_height_m()` (`Σ h_i`, ground → lid).

`auto_dx` itself is untouched — `auto_dx_stackup` is a new function, not a
refactor of the existing one, so no behavior change reaches any existing caller.

Unit tests (`crates/yee-engine/src/automesh.rs`'s `mod tests`): each rule binding
in turn (wavelength via a higher-ε_r second layer; h/3 via a thin layer among a
thick one; feature/2 via a narrow-gap fixture; lid b/16 via
`Stackup::symmetric_stripline(4.4, 3.2e-3)` — the same numbers as the stripline
gates — under a loose wavelength/h-3/feature budget, asserting the other three
terms are provably looser, not just that the number matches); a clamp case
(mirrors `auto_dx_is_clamped`); and a single-layer/no-lid `Stackup` built from
`layout.substrate`'s own fields reproducing `auto_dx`'s result **bit-for-bit**
(`assert_eq!`, not a tolerance) across three `f_max_hz` values — the spec's
consistency check.

### 2. Gate `engine-automesh-stackup-001` (Task 2, commit `fed0a87`)

New file `crates/yee-engine/tests/automesh_stackup.rs`. Rebuilds the
`engine-stripline-eeff-001` fixture (identical `EPS_R = 4.4`, `B_M = 3.2e-3`,
`W_M = 1.5e-3`, `F0_HZ = 5.0e9`, `PORT_R_OHM = 50.0`, same trace/port/probe layout)
with **no hand-set dx anywhere**:

- `dx_seed = auto_dx_stackup(&layout, &stack, f_max_hz)` seeds the grid
  (`f_max_hz = F0_HZ + bw/2` from the drive's own Gaussian half-bandwidth — the
  "f_max from the drive" idiom `board_automesh.rs` already uses).
- The CPML margin is held as a **metres** constant (`MARGIN_M = 4.0e-3`, the same
  physical 4 mm the hand-tuned gate picked as 20 cells at its hand-set 0.2 mm dx)
  and converted to cells via `(MARGIN_M / dx_seed).round()` — cell-denominated
  quantities derive from the dx the rulebook actually returns, the ADR-0204
  constant-physics hygiene applied to a second gate family.
- The test recomputes all four rule terms itself (wavelength λ/20, per-layer h/3,
  feature/2, lid b/16) rather than hardcoding an expectation, takes the min,
  `eprintln`s which term binds plus every term, and `assert_eq!`s the binding
  rule's name is `"lid b/16 (ADR-0215)"` (the assert message names which rule
  bound instead, and why a different one binding would matter, if it doesn't).
  Separately asserts the returned `dx_seed` matches the recomputed lid term to
  `< 1e-12`.
- The rest — `voxelize_stackup` → `JobSpec` → submit → time-gated two-probe
  single-bin DFT phase advance → v_p → ε_eff — is unchanged from
  `stripline_eeff.rs`, reading `dx`/`margin_cells` instead of the old
  `DX_M`/`MARGIN_CELLS` constants. Tolerance: the same `rel_err ≤ 0.02` bar as
  `engine-stripline-eeff-001` (not widened).

`stripline_eeff.rs` itself is untouched.

## Measured result

```
engine-automesh-stackup-001: dx = 0.2000 mm, binding rule = lid b/16 (ADR-0215)
  (0.2000 mm) [wavelength 1.1910, h/3 0.5333, feature/2 0.7500, lid b/16 0.2000]
  grid 1184x48x16, trace at k = 8 (b = 16 cells), L = 228.7 mm
  eps_eff = 4.4029 vs exact TEM 4.4000 -> err 0.065 % (dphi = 2.1109 rad over
  9.600 mm, 7027 steps, gate 6827)
test automesh_stackup_matches_the_exact_tem_value ... ok
finished in 22.64s – 22.92s
```

The lid rule binds exactly as the spec predicted (0.2000 mm, the loosest of the
other three terms by ≥ 2.65×), and the rulebook-seeded run reproduces the
hand-tuned `engine-stripline-eeff-001` gate's dx (0.200 mm), grid shape
(1184×48×16), and ε_eff error (**0.065 %**) — bit-for-bit, not merely "close" —
because both grids are the same physical dx applied to the same fixture.

## Tolerance pinned

`rel_err ≤ 0.02` (2 %), matching `engine-stripline-eeff-001`'s bar unchanged.
Measured 0.065 % gives ~30× headroom, identical to the hand-tuned gate's own
margin — no widening, no narrowing.

## Bit-exactness / regression discipline (unmodified gates, every commit)

The binding gate command — `cargo test -p yee-compute --release --test
graded_uniform_bitexact --test gpu_graded_parity --test gpu_cpu_parity --
--include-ignored` — stayed green (5/5) after both Task 1 and Task 2 commits, GPU
evidence confirmed real (`compute-002: running on adapter 'NVIDIA GeForce RTX 5060
Ti'`, not SKIPPED; family-rel L2/L∞ ~1e-7–6e-7 both runs). `cargo test -p
yee-engine --release` (default set) reported **45/45 unit tests pass, 0 failed**
every run after both commits; a pre-existing flaky GPU/wgpu-teardown SIGSEGV
*after* `test result: ok` (reproduced on the pre-Task-1 commit too, gone under
`--test-threads=1`) is unrelated to this track's diff and out of lane — recorded
in the Task 1/2 reports, not fixed here. The three stripline gates re-ran green
and unmodified after both commits (`git diff` on `stripline_eeff.rs` /
`stripline_z0.rs` / `stripline_alpha.rs` empty throughout): `stripline_alpha`
72–83 s, `stripline_eeff` 23–24 s, `stripline_z0` ×2 30–38 s, ~125–145 s total —
within the ~3 min budget. Workspace clippy (default + `--no-default-features` on
`yee-compute`) and `cargo fmt --check --all` clean before every commit;
`missing_docs` clean.

## Verdict

**GO.** `auto_dx_stackup`'s rulebook alone — with no fixture-specific tuning —
reproduces the hand-picked `engine-stripline-eeff-001` dx bit-for-bit and lands
inside the same certified ≤ 2 % tolerance on the first honest run. The ADR-0215
b/16 lesson is now load-bearing code, not just a prose warning: the lid rule is
the term that actually binds on the reference fixture, exactly as predicted.
Push-button meshing (FS.0's walking skeleton) now covers stackup boards.

## What remains of FS.4.2 (queued, not attempted here)

Per the spec's non-goals: a graded `auto_spacings` stackup variant (uniform-dx
rulebook first, per "walking skeleton first" — graded multilayer is a follow-on);
a MoM cross-check of the stackup/dielectric-loss models; automesh awareness of
vias or arcs in a stackup. None of these were scoped into FS.4.2c.
