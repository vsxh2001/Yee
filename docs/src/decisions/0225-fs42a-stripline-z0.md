# ADR-0225: FS.4.2a — stripline Z₀: H-field probes + V/I extraction gate

**Date:** 2026-07-23 · **Status:** accepted · **Track:** FS.4 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-23-fs4-2a-stripline-z0-design.md`
**Plan:** `docs/superpowers/plans/2026-07-23-fs4-2a-stripline-z0.md`
**Predecessor:** FS.4.0 (ADR-0215) — `voxelize_stackup` + `engine-stripline-eeff-001` (ε_eff
0.065 % vs exact TEM); FS.4.1 (ADR-0221) — through/blind vias, `engine-stackup-via-001`.

## Context

FS.4.0/4.1 validated stripline **propagation** (ε_eff, notch frequencies). The FS.4 roadmap row's
next named gate is stripline **Z₀** — the first *impedance* (not propagation-constant) full-wave
validation in the engine flow. Z₀ = V/I needs a measured voltage and a measured current on the
line. V is a column of `Ez` probes summed × Δz — already possible. I is an Ampère loop
`∮H·dl` around the trace — the engine had no H-field probes (`Probe.component` is
`EComponent`-only, E-field probes only). That gap is Task 1's deliverable; Task 2 is the gate
itself.

## Decision

### 1. H-field probes (Task 1, commit `758d206`)

`Drive` gains a parallel `h_probes: Vec<HProbe>` field (`HProbe { component: HComponent, cell }`,
`HComponent::{Hx,Hy,Hz}`) — **not** a widened `Probe.component`/`EComponent`. Chosen over widening
because several existing tests construct `Drive`/`Probe` via full struct literals; widening
`EComponent` would have forced edits to `graded_uniform_bitexact.rs` and `gpu_graded_parity.rs`,
two of the three files this track's binding verify command requires to stay byte-**unmodified**.
The parallel-field shape leaves `Probe`/`EComponent` completely untouched — zero edits to any
pinned gate file. CPU recording is exact (`h_probe_series` on `CpuFdtd`, mirroring
`probe_series`); GPU recording was implemented too (a contained, append-only addition to the
`record_probes` WGSL pipeline and `drv_idx`/`drv_data` layout — H is already packed in the same
field arena right after E) rather than rejected as `Unsupported`, since the edit stayed small.

**Timing is load-bearing and is documented on `HProbe` itself**: the H state recorded alongside a
given E-probe sample is a *half step behind* it — `update_h` runs at the top of each iteration
(t = (n+½)·Δt) before the co-recorded E sample is written at t = (n+1)·Δt. Task 2's DFT phasors
use each series' own true sample time rather than pairing same-index samples as simultaneous.

**Regression discipline**: 7 existing test files needed a one-line `h_probes: vec![]` addition to
their full-literal `Drive` constructions (behavior-neutral, empty vec). The existing
`cpu_drive_parity.rs::driven_step_is_bit_exact_against_reference` — otherwise untouched — still
asserts the E-probe series matches the `yee_fdtd` reference bit-for-bit, confirming E-probe
behavior is provably unchanged.

**GPU parity** (`gpu_h_probe_parity.rs`, real hardware, `NVIDIA GeForce RTX 5060 Ti`): H-probe
rel L2 = 9.5e-7 / 1.0e-6 on two probes, rel L∞ ≤ 6.0e-7 — same FP32 idiom as `gpu_cpu_parity`
(compute-002).

### 2. Gate `engine-stripline-z0-001` (Task 2, commit `1ea414a`)

New file `crates/yee-engine/tests/stripline_z0.rs`. Constructs
`yee_compute::{FdtdSpec, Fields, Materials, Boundary, Drive, CpuFdtd}` directly rather than going
through `yee_engine::JobSpec::submit` — `JobSpec`/`ProbeSpec` carry E probes only, and every
`JobSpec` construction workspace-wide (29 call sites across 4 crates) is a full struct literal
with no `Default`, so widening it for a single gate's need would force unrelated edits across all
of them. This uses the same primitives `run_job` itself calls internally, from `yee-engine`'s own
(already-a-dependency) test suite — in-lane, zero edits to `yee-engine/src/`.

**Fixture**: symmetric stripline, ε_r = 2.2, b = 3.2 mm (16 cells at dx = 0.2 mm — the
ADR-0215/0221 confined-mode lesson), w = 2.6 mm (13 cells, w/b = 0.8125) — solved by bisection for
Z₀_exact = 50 Ω and rounded to the nearest whole cell, not guessed. Box margin 20 cells (TE₁₀
cutoff ≈ 9.5 GHz, clear of the ≈ 8.4 GHz drive-band top). L = 8·λ_g. Grid 1388×53×16, 6035 steps,
~24 s release.

**Closed form**: `Z₀ = (η₀/4√ε_r)·K(k′)/K(k)`, `k = tanh(πw/2b)`, `k′ = sech(πw/2b)`, `K` via a
~10-line AGM iteration (Abramowitz & Stegun 17.6). Cross-checked in an always-on test against the
independent Wheeler/Cohn fit (`30π/√ε_r/(w/b+0.441)`, Pozar §3.8): **0.052 % agreement** at the
fixture's w/b, and the physical-sanity check that Z₀ falls as w/b grows.

**A spec bug caught before running any FDTD**: the design spec's own text labels
`k = sech(πw/2b), k′ = tanh(πw/2b)` — swapped from the standard/Wikipedia "Stripline" convention.
Plugging the spec's literal labels into `Z₀ = K(k′)/K(k)` makes Z₀ *increase* with w/b, which
contradicts basic transmission-line physics (wider trace ⇒ more capacitance per length ⇒ *lower*
Z₀, Z₀ → ∞ as w → 0) and disagreed with the Wheeler fit by tens of percent instead of the fit's
usual ≲ 1–2 %. Re-deriving with the standard labelling (`k = tanh`, `k′ = sech`, same `K(k′)/K(k)`
formula shape) restores the correct decreasing trend and the < 0.1 % Wheeler agreement. The
shipped code uses the corrected convention, documented inline in the module doc's "Closed form"
section and via the always-on cross-check test — this is the same "verify the reference, don't
blindly trust an unverified formula" discipline the track's binding constraints require, applied
one level up: to the design spec's own text, before it could corrupt the gate.

**V(t)**: column of `Ez` probes from the ground plane (k=0) up to (excluding) the trace plane
(k=k_trace), summed × Δz, at a plane 2.5 guided wavelengths downstream of the port — past the
launch transient, same plane-placement hygiene as `engine-stripline-eeff-001`. Symmetric-stripline
argument for why the lower-half integral equals the full line voltage is in the module doc.

**I(t)**: a rectangular Ampère loop, one Δz tall, straddling the trace plane exactly (`Hy` at
`k_trace−1` below, `Hy` at `k_trace` above), spanning the trace's y-extent plus a 5-cell guard
(`Hz` at the two side legs). Derived from first principles (not assumed) that this loop *is* the
FDTD curl loop around `Ex(i, j, k_trace)` — the Ampère-law surface it bounds is pierced only by
the trace's own PEC surface current — and that summing adjacent unit loops telescopes the interior
`Hz` terms away, leaving

```
I = Δy·[Σ_j Hy(j,k_trace−1) − Σ_j Hy(j,k_trace)] + Δz·[Hz(j_hi,k_trace) − Hz(j_lo−1,k_trace)]
```

Because Ampère's law is exact and total current is conserved (∇·(J+∂D/∂t)=0 identically), the
guard-cell choice is provably not load-bearing for correctness — any non-crossing loop around the
trace measures the same I.

**Staggering, handled explicitly (not ignored)**:

- **Time**: Task 1's documented half-step offset — `Ez` at t=(m+1)·Δt, `Hy`/`Hz` at t=(m+½)·Δt.
  The gated single-bin DFT phase-references each series at its own true sample time, so this is
  handled exactly, no interpolation.
- **Space**: `Ez`'s probe index is an integer-x node; `Hy`/`Hz`'s same numeric index is the
  half-x node a half cell downstream. Not corrected — quantified instead: at this fixture's f0/λ_g
  the induced phase error is β·Δx/2 ≈ 0.019 rad (≈ 1.1°), a cos magnitude error < 0.02 % — three
  orders under the gate's tolerance.

**Z₀ = |V(f0)|/|I(f0)|** from the gated phasor magnitudes (a time-domain single-bin DFT, same
idiom as `engine-stripline-eeff-001`'s phase extraction).

## Measured result

```
$ cargo test -p yee-engine --release --test stripline_z0 -- --ignored --nocapture
engine-stripline-z0-001: grid 1388x53x16, trace at k = 8 (b = 16 cells), w/b = 0.8125,
  Z0_exact = 50.6651 Ohm, L = 269.5 mm
  Z0_meas = 50.0209 Ohm vs exact 50.6651 Ohm -> err 1.271 % (|V| = 2.333e0, |I| = 4.664e-2,
  6035 steps, gate 5835)
test stripline_z0_matches_the_exact_closed_form ... ok   (24.24 s)
```

**1.271 % measured, first run, no root-cause detour needed** — well inside the plan's ≤ 5 % target
and never approached the 10 % STOP-and-root-cause threshold. `V` and `I` are guarded against a
silent-zero wiring bug with a 1e-3 non-triviality floor, three orders under the smaller measured
magnitude (|I| ≈ 0.047).

## Tolerance pinned

`rel_err ≤ 0.05` (5 %, matches the plan's explicit target). Measured 1.271 % gives ~4× headroom.
Never widened — not needed, since the first measurement already cleared tolerance.

## Bit-exactness / regression discipline (unmodified gates, every commit)

The binding gate command —
`cargo test -p yee-compute --release --test graded_uniform_bitexact --test gpu_graded_parity
--test gpu_cpu_parity -- --include-ignored` — stayed green after both commits, with all three
gate files confirmed byte-unmodified (`git status --porcelain` empty on each). `engine-stripline-eeff-001`
re-run after Task 2 and unchanged (ε_eff 0.065 %, the pre-existing pinned value). Both crates'
workspace clippy (default + `--no-default-features` on `yee-compute`) and `cargo fmt --check --all`
clean before every commit; `missing_docs` clean.

## Verdict

**GO.** The H-probe machinery is a clean, non-invasive addition (parallel field, zero churn to
pinned gates); the Z₀ gate passed its first honest measurement at 1.271 %, confirming the
extraction method (V-column integral + Ampère-loop current, staggering handled per-series) is
physically sound on a symmetric stripline. This is the first engine-flow gate that measures an
*impedance* rather than a propagation constant — the V/I plumbing (H-probes, gated-DFT phasor
extraction) is now reusable for any future port-impedance or de-embedding work.

## What remains of FS.4.2 (queued, not attempted here)

Per the design spec's non-goals: per-layer tan δ in the engine materials model, MoM multilayer
cross-check, automesh awareness of stackups (the FS.4.0-flagged "≥ 16 cells across b" rule folded
into the rulebook), microstrip Z₀ (harder reference — quasi-TEM, not exact closed-form, and the
mom-002 port ADR-0064 lesson means a numerical-port frame mapping is not a shortcut here either),
S-parameter port-impedance renormalization to a measured Z₀ (would consume `yee_io::network`'s
FS.6.1 `renormalize`, once a non-50 Ω measured Z₀ shows up in a real fixture).

If H-probes are ever needed through the `yee_engine::JobSpec` job protocol (studio/server
exposure rather than a test-only direct `yee-compute` construction): `ProbeSpec.component` would
need widening to accept `"hx"/"hy"/"hz"` plus a parallel `JobResult.h_probes` populated from both
run paths — not attempted here, out of scope for a single gate.
