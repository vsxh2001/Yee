# Tutorial 17 — Designing an RF filter

This tutorial walks the first stages of Yee's end-to-end filter-design
flow: you write a typed **filter specification**, synthesize a lowpass
prototype + coupling matrix from it, sweep the ideal response, and grade
that response against a **spec mask** — all without touching an EM
solver. You will do it twice: once from the `yee filter synth` CLI
(scriptable, CI-friendly), and once in the `yee-studio` desktop app
(interactive, live-updating). By the end you will have a Touchstone
`.s2p` you can open in any S-parameter tool and a PASS/FAIL verdict
against your requirements.

## The design flow

Yee's stated final goal is an interactive RF-filter designer shipped as
**both a desktop app and a web app** — one `egui`/`eframe` codebase that
targets native and WASM (see `FILTER-DESIGN-ROADMAP.md` and ADR-0089).
The full pipeline runs spec → prototype → coupling matrix → physical
dimensions → parametric layout → full-wave-verified S-parameters →
fabrication files. The flow is *synthesis-assisted interactive*: the
tool proposes each stage, you inspect and approve, then it proceeds.

What is shipped **today** — and what this tutorial covers — is the
"light flow" front-end:

- **`yee-synth`** — closed-form prototypes (Butterworth, Chebyshev),
  the bandpass transform, and the all-pole coupling matrix + external Q.
- **`yee-filter`** — the `FilterSpec` data model, the closed-form ideal
  response, and the spec-mask gate (`check_mask`).
- **`yee filter synth [--plot]`** — the CLI front-end over those crates.
- **`yee-studio`** — the `eframe` desktop app that wraps the same logic
  in live spec-editor / synthesis / plot panels.

The heavy stages — **EM-in-the-loop dimensional synthesis** (mapping the
abstract coupling matrix to physical resonator lengths and gaps via
coupling extraction + `yee-surrogate` BO with FDTD in the loop, Phase
F1.1+), layout generation, full-wave verification, and fabrication
export — are forthcoming. This tutorial is the part you can run end to
end on a laptop with nothing but the workspace checkout.

## Goal

Take a 0.5 dB-ripple Chebyshev bandpass spec (2 GHz centre, 10 %
fractional bandwidth, order 5, with a 40 dB rejection point at 2.4 GHz),
synthesize it, confirm the spec-mask verdict is **PASS**, and write a
Touchstone `.s2p`. Then open the same design in `yee-studio` and watch
the response and verdict update live as you drag the order and mask
sliders.

## Prerequisites

- Rust 1.92+ and a workspace checkout. No Gmsh, no Python, no CUDA.
- For the `--plot` flag and the Studio app: the plotters / egui native
  deps. On Linux that means `libfontconfig1-dev` and `pkg-config`
  (`sudo apt install libfontconfig1-dev pkg-config`). The bare
  `yee filter synth` (no plot) needs neither.

## Write a FilterSpec

A `FilterSpec` is the design intent: the response class, the
approximation (response shape), the centre frequency and fractional
bandwidth, an optional explicit order, the system impedance, and the
spec **mask** the response is graded against. It is plain TOML. The
repo ships a satisfiable example at
`crates/yee-cli/tests/fixtures/cheb_bpf.toml`:

```toml
response = "Bandpass"
f0_hz = 2.0e9
fbw = 0.10
order = 5
z0_ohm = 50.0

[approximation.Chebyshev]
ripple_db = 0.5

[mask]
passband_ripple_db = 0.5
return_loss_db = 9.0
# Each stopband row is [frequency_hz, minimum_rejection_db].
stopband = [[2.4e9, 40.0]]
```

Field by field:

- **`response`** — the frequency-response class. Today's synthesis
  evaluates `"Bandpass"`; `"Lowpass"`, `"Highpass"`, and `"Bandstop"`
  are reserved for later phases.
- **`approximation`** — the response shape. `[approximation.Chebyshev]`
  with a `ripple_db` field gives an equiripple passband;
  `approximation = "Butterworth"` (a bare string, no fields) gives a
  maximally-flat one.
- **`f0_hz`** — centre frequency, Hz.
- **`fbw`** — fractional bandwidth `(f2 − f1) / f0`. Here `0.10` puts
  the band edges at 1.9 / 2.1 GHz.
- **`order`** — explicit filter order N. Omit it (`order` absent) to let
  Yee estimate the minimum order that meets the worst-case stopband
  point.
- **`z0_ohm`** — system reference impedance, written into the Touchstone
  option line.
- **`[mask]`** — `passband_ripple_db` (max allowed in-band insertion-loss
  ripple), `return_loss_db` (min required in-band return loss; a *larger*
  value is stricter), and `stopband` (a list of `[frequency_hz,
  min_rejection_db]` rows, each requiring `|S21|` at that frequency to be
  at least that many dB down).

> **A note on the return-loss value.** A 0.5 dB-ripple Chebyshev caps
> in-band return loss at ≈ 9.64 dB, so the mask asks for `9.0` dB, not
> `10.0` — a stricter `return_loss_db` than the shape can deliver would
> make the mask unsatisfiable for this filter. Tighten the ripple (or
> raise the order) if you need a better match.

## Synthesize from the CLI

Point `yee filter synth` at the spec:

```bash
cargo build --release -p yee-cli
./target/release/yee filter synth crates/yee-cli/tests/fixtures/cheb_bpf.toml
```

The handler parses the spec, synthesizes the prototype and coupling
matrix, sweeps the closed-form ideal response over a 401-point grid
spanning `f0·(1 ± 6·fbw/2)`, grades it against the mask, and writes a
Touchstone `.s2p` next to the spec (here
`crates/yee-cli/tests/fixtures/cheb_bpf.s2p`). You will see:

```
Filter synthesis (Chebyshev { ripple_db: 0.5 }, order N=5)
  f0 = 2.000000e9 Hz   FBW = 0.1000   Z0 = 50 Ohm
  prototype g-values: g0=1.0000  g1=1.7058  g2=1.2296  g3=2.5409  g4=1.2296  g5=1.7058  g6=1.0000
  external Q: Qe_in=17.0582  Qe_out=17.0582
  coupling matrix M (normalized, 5x5):
    [ +0.0000  +0.6905  +0.0000  +0.0000  +0.0000 ]
    [ +0.6905  +0.0000  +0.5657  +0.0000  +0.0000 ]
    [ +0.0000  +0.5657  +0.0000  +0.5657  +0.0000 ]
    [ +0.0000  +0.0000  +0.5657  +0.0000  +0.6905 ]
    [ +0.0000  +0.0000  +0.0000  +0.6905  +0.0000 ]
  mask: passband ripple 0.499 dB (spec 0.500), in-band RL 9.641 dB (spec 9.000)
  stopband 2.4000e9 Hz: rejection 70.54 dB (required 40.00 dB) OK
  wrote Touchstone: crates/yee-cli/tests/fixtures/cheb_bpf.s2p
VERDICT: PASS
```

The process exits **0** on PASS and **1** on FAIL, so the command drops
straight into a CI gate or a shell `&&` chain.

Reading the report top to bottom:

- **prototype g-values** — the order-5 lowpass-prototype element values
  (`g0 … g6`). These are the published Chebyshev 0.5 dB N=5 values; they
  are the seed every downstream stage maps from.
- **external Q** — `Qe_in` / `Qe_out`, the input/output coupling that
  ties the end resonators to the 50 Ω terminations.
- **coupling matrix M** — the normalized symmetric tridiagonal
  coupled-resonator matrix (zero diagonal ⇒ synchronously tuned). The
  off-diagonal entries are the inter-resonator couplings.
- **mask** — the achieved in-band ripple and return loss versus the spec
  values, then one line per stopband point with the achieved rejection,
  the requirement, and `OK` / `UNDER`.
- **VERDICT** — `PASS` iff in-band ripple ≤ spec, in-band RL ≥ spec, and
  every stopband point meets its rejection.

### Render the response with the spec mask

Add `--plot` to also draw `|S21|` (dB) with the mask's forbidden regions
shaded — a passband floor at `−passband_ripple_db` over the band edges,
and a ceiling at `−rejection` over a ±2 % band around each stopband
point:

```bash
./target/release/yee filter synth \
    crates/yee-cli/tests/fixtures/cheb_bpf.toml \
    --output cheb_bpf.s2p --plot cheb_bpf.png
```

`--output` overrides the Touchstone path; `--plot` takes a `.png` or
`.svg` target (the extension picks the format). The extra stdout line is
`wrote plot: cheb_bpf.png`. The plot makes the verdict visual: a passing
design keeps its `|S21|` trace out of every shaded box.

If you want to drive the same synthesis as a registered validation gate
rather than a one-off, `yee validate synth` runs the Filter Phase F0
gates (`synth-*` / `filt-*`).

## The Studio app

`yee-studio` is the interactive front-end over exactly the same
`yee-filter` logic. Launch the native desktop build:

```bash
cargo run -p yee-studio
```

The crate's default `desktop` feature pulls in the `eframe`/`egui`/`wgpu`
windowing shell; the app opens seeded with the same satisfiable
Chebyshev 0.5 dB N=5 bandpass design used above. Three regions track the
spec live:

- a left **spec-editor** panel — drag `f0`, FBW, order, ripple, return
  loss, and stopband points; pick Chebyshev or Butterworth;
- a central **synthesis** panel — the prototype g-values, the coupling
  matrix grid, the external Q, and a coloured **PASS** (green) / **FAIL**
  (red) verdict with the same per-line mask notes the CLI prints;
- an `egui_plot` **`|S21|`-vs-mask** view — the response trace with each
  forbidden mask region shaded on its violating side.

Every edit re-runs `synthesize` → `ideal_response` → `check_mask`, so
you can watch a design slide from PASS to FAIL as you drop the order
(try N=2 against the 40 dB stopband) and back as you raise it.

Under the hood the design state lives in a `StudioState` value that is
deliberately **egui-free and WASM-safe** (ADR-0090 / ADR-0092): only
`app.rs` and the binary entry depend on `egui`/`eframe`, gated behind the
`desktop` Cargo feature. A `--no-default-features` build compiles just
the flow logic with no GUI in the dependency graph — the groundwork for
the forthcoming browser build (App.1), where this same light flow runs
fully client-side in WASM with no server.

## What you get, and what's next

From this stage you have two concrete artifacts:

- a Touchstone **`.s2p`** of the ideal response — open it in any
  S-parameter tool (or feed it straight back through Yee's `yee-io`
  Touchstone reader);
- a spec-mask **PASS/FAIL** verdict — a machine-checkable gate on whether
  the synthesized response meets your requirements.

Note that the `.s2p` here is the *ideal closed-form* response — a
magnitude model from the prototype, not yet an EM result. The next
stages of the flow turn that abstract design into a real, manufacturable
filter: **dimensional synthesis** maps the coupling matrix to physical
resonator dimensions (coupling extraction + `yee-surrogate` Bayesian
optimization with the FDTD solver in the loop), **layout generation**
produces the parametric geometry, **full-wave verification** EM-simulates
the complete layout and re-checks it against this same spec mask, and
**fabrication export** writes KiCad/Gerber (planar/lumped) or STEP
(waveguide). Those stages are tracked in `FILTER-DESIGN-ROADMAP.md`
(Phase F1.1+ and the App/Studio track); this tutorial is the front of
that pipe.

## Next

For the EM solvers the verification stage will lean on, see
[Tutorial 3 — FDTD cavity resonance](03-fdtd-cavity.md) and
[Tutorial 11 — FDTD cavity resonance from
Python](11-fdtd-cavity-resonance-from-python.md). For the
natural-language route into a spec, see [Tutorial — Natural-language
design surface](04-nl-design-surface.md).
