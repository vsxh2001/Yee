# Tutorial 17 — Designing an RF filter

This tutorial walks Yee's end-to-end filter-design flow: you write a typed
**filter specification**, let the tool **recommend a realization technique**,
synthesize a coupling matrix / ladder, **compare techniques**, pick one of six
shipped microstrip / lumped realizations, get **physical dimensions** and a
**parametric board layout**, sanity-check **tolerance / yield** and a **bill of
materials**, and **export Gerber / KiCad** fabrication files. You will do it two
ways: in the interactive **`yee-studio-web`** studio (a pure-Rust Dioxus app that
runs in the browser), and from the scriptable **`yee filter synth`** CLI.

One thing is honest up front: the studio's **Verify** stage grades at the
**circuit / synthesis level** (the synthesized ideal — or, for the lumped flow,
the *realized* LC ladder — response against your spec mask). **Full-wave EM
verification of the finished board is not run here** — it is a separate native
step and a deferred research frontier (ADR-0133). The tool never hides that gap;
neither does this tutorial.

## The design flow

Yee's filter designer is an interactive, synthesis-assisted flow: the tool
proposes each stage, you inspect and refine, then it proceeds. The full pipeline
is:

```
Spec → (recommender suggests a technique) → Synthesis (coupling matrix / ladder)
     → Compare techniques → pick a realization → Dimensions + parametric layout
     → Tolerance / yield + BOM (lumped) → Gerber / KiCad export
     → Verify (circuit-level vs the spec mask)
```

Everything in this tutorial runs **on a laptop with nothing but the workspace
checkout** — no Gmsh, no CUDA, no Python. The whole design flow is intentionally
serde-only and WASM-safe (ADR-0089/0099), which is exactly why the same code
ships as a browser app.

### The two front-ends

- **`yee-studio-web`** — the interactive studio. A pure-Rust [Dioxus] app that
  compiles to WebAssembly and runs **client-side in the browser**. It is
  deployed live on GitHub Pages:

  **<https://vsxh2001.github.io/Yee/studio/>**

  This is *the* studio. (The earlier `eframe`/`egui` desktop app, `yee-studio`,
  was **retired** in App.D.2 / ADR-0130. The `yee-studio` crate still exists, but
  only as a headless `StudioState` logic library — there is no longer a desktop
  app binary to run. Do not look for one.)

- **`yee filter synth`** — the CLI front-end over the same `yee-filter` engine.
  Scriptable and CI-friendly: it synthesizes, grades against the spec mask,
  writes a Touchstone `.s2p`, dimensions an edge-coupled board, and can emit a
  layout SVG, Gerber, and KiCad PCB.

Both drive the same crates:

- **`yee-synth`** — closed-form prototypes (Butterworth, Chebyshev), the
  bandpass transform, the all-pole coupling matrix + external Q.
- **`yee-filter`** — the `FilterSpec` data model, the ideal response, the
  spec-mask gate (`check_mask`), the technique **recommender**
  (`recommend_technique`), the six dimensioners (`dimension_edge_coupled`,
  `dimension_hairpin`, `dimension_combline`, `dimension_interdigital`,
  `dimension_stepped_impedance`, plus the lumped `synthesize_lumped`), and the
  tolerance / BOM helpers (`monte_carlo_yield`, `select_components`).
- **`yee-export`** — the Gerber (`F.Cu` + `Edge.Cuts`) and KiCad (`.kicad_pcb`)
  board emitters.

[Dioxus]: https://dioxuslabs.com/

## Goal

Take a 0.5 dB-ripple Chebyshev bandpass spec (2 GHz centre, 10 % fractional
bandwidth, order 5, with a 40 dB rejection point at 2.4 GHz), recommend a
technique, synthesize it, confirm the spec-mask verdict is **PASS**, write a
Touchstone `.s2p`, get edge-coupled microstrip dimensions, and export a Gerber +
KiCad board — all from the CLI. Then open the same design in the studio and
watch it re-derive live as you edit the spec, switch techniques, and step
through the stages.

## Prerequisites

- Rust 1.92+ and a workspace checkout. No Gmsh, no Python, no CUDA.
- For the CLI `--plot` flag: the plotters native deps. On Linux that means
  `libfontconfig1-dev` and `pkg-config`
  (`sudo apt install libfontconfig1-dev pkg-config`). The rest of
  `yee filter synth` (no `--plot`) needs neither.
- To run the **studio locally** (optional — it is already live on the URL above):
  the Dioxus CLI and the wasm target:

  ```bash
  rustup target add wasm32-unknown-unknown
  cargo install dioxus-cli --version 0.6.3 --locked
  cd crates/yee-studio-web && dx serve --platform web
  ```

## Write a FilterSpec

A `FilterSpec` is the design intent: the response class, the approximation
(response shape), the centre frequency and fractional bandwidth, an optional
explicit order, the system impedance, and the spec **mask** the response is
graded against. It is plain TOML. The repo ships a satisfiable example at
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

- **`response`** — the frequency-response class: `"Bandpass"`, `"Lowpass"`,
  `"Highpass"`, or `"Bandstop"`. The recommender reads all four; the live
  synthesis flows build band-pass (every distributed technique except
  stepped-impedance, plus lumped) and low-pass (stepped-impedance).
- **`approximation`** — the response shape. `[approximation.Chebyshev]` with a
  `ripple_db` field gives an equiripple passband; `approximation =
  "Butterworth"` (a bare string, no fields) gives a maximally-flat one.
- **`f0_hz`** — centre frequency (band filters) or cutoff (low/high-pass), Hz.
- **`fbw`** — fractional bandwidth `(f2 − f1) / f0`. Here `0.10` puts the band
  edges at 1.9 / 2.1 GHz.
- **`order`** — explicit filter order N. Omit it to let Yee estimate the minimum
  order that meets the worst-case stopband point.
- **`z0_ohm`** — system reference impedance, written into the Touchstone option
  line.
- **`[mask]`** — `passband_ripple_db` (max allowed in-band insertion-loss
  ripple), `return_loss_db` (min required in-band return loss; a *larger* value
  is stricter), and `stopband` (a list of `[frequency_hz, min_rejection_db]`
  rows, each requiring `|S21|` at that frequency to be at least that many dB
  down).

> **A note on the return-loss value.** A 0.5 dB-ripple Chebyshev caps in-band
> return loss at ≈ 9.64 dB, so the mask asks for `9.0` dB, not `10.0` — a
> stricter `return_loss_db` than the shape can deliver would make the mask
> unsatisfiable for this filter. Tighten the ripple (or raise the order) if you
> need a better match.

## Let the tool recommend a technique

Before committing to a topology, `yee-filter` can suggest one. The recommender
(`recommend_technique`, ADR-0136) is a **deterministic decision tree** keyed on
the response class, the centre / cutoff frequency, and the fractional bandwidth
— with the thresholds (a ~500 MHz distributed-feasibility floor, the 5 % / 20 %
fractional-bandwidth bands) pinned by the `tech_001` gate so they cannot drift.
It returns a primary technique, a plain-language **rationale** that names the
deciding factor, and a ranked list of **alternatives**, each with a one-line
tradeoff note.

For the spec above (band-pass, 2 GHz, 10 % FBW) the tree lands on
**edge-coupled** — distributed, ≥ 500 MHz, moderate 5–20 % bandwidth — and
offers **hairpin** as the compact alternative. The shape of the tree:

- **Low-pass:** cutoff ≥ 500 MHz → stepped-impedance; else lumped LC.
- **High-pass:** lumped LC (Yee's distributed techniques are LP/BP-oriented —
  the recommender says so honestly).
- **Band-pass / band-stop:** below 500 MHz → lumped LC; ≥ 500 MHz with FBW ≥ 20 %
  → edge-coupled; 5–20 % → edge-coupled (hairpin alternative); < 5 % →
  interdigital (combline alternative).

In the studio this drives the guided entry point: type a frequency, bandwidth,
and response class, and the **Technique** stage shows the recommendation plus the
ranked alternatives, each as a clickable card that routes into its flow.

## The six realizations

The studio's coupled-resonator gallery is **complete** — all six techniques are
live and routable (no "coming soon" placeholders remain). Each maps to a real
dimensioner in `yee-filter` and produces a microstrip / lumped board you can
export:

| Technique | Response | Engine | Notes |
|---|---|---|---|
| **Edge-coupled** | Band-pass | `dimension_edge_coupled` | Parallel-coupled half-wave lines; the broad default. |
| **Hairpin** | Band-pass | `dimension_hairpin` | U-folded half-wave resonators — same synthesis, smaller board. |
| **Combline** | Band-pass | `dimension_combline` | Grounded quarter-wave resonators with end-loading capacitors. |
| **Interdigital** | Band-pass | `dimension_interdigital` | λg/4 lines short-circuited at *alternating* ends (no loading cap); compact, high-Q. |
| **Lumped LC** | Band-pass | `synthesize_lumped` | Discrete L/C ladder → SMD parts + BOM + tolerance/yield. |
| **Stepped-impedance** | Low-pass | `dimension_stepped_impedance` | Alternating high-/low-impedance line sections. |

The studio's **left rail** adapts to the chosen technique. The five distributed
band-pass / low-pass techniques run a **six-stage** rail — Spec → Technique →
Synthesis → Layout → Verify → Export. The **lumped** flow runs an **eight-stage**
rail that inserts **Components** (BOM) and **Tolerance** (yield) between Synthesis
and Layout.

Before picking, the **Compare** view (`compare_techniques`, ADR-0142) tabulates
each technique's graded mask metrics side by side, and the **response overlay**
(`overlay_curves`, ADR-0143) plots their `|S21|` traces against the spec mask on
one axis — so the technique choice is a measured comparison, not a guess.

## Synthesize from the CLI

The CLI is the fastest way to drive the whole flow non-interactively. Point
`yee filter synth` at the spec:

```bash
cargo build --release -p yee-cli
./target/release/yee filter synth crates/yee-cli/tests/fixtures/cheb_bpf.toml
```

The handler parses the spec, synthesizes the prototype + coupling matrix, sweeps
the closed-form ideal response over a 401-point grid spanning `f0·(1 ± 3·fbw)`,
grades it against the mask, writes a Touchstone `.s2p` next to the spec, and then
dimensions an **edge-coupled** microstrip board on the substrate
(FR-4 defaults `εr = 4.4`, `h = 1.6 mm`; override with `--eps-r` / `--h-mm`).
You will see:

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
  substrate: eps_r = 4.4000   h = 1.6000 mm
  physical dimensions (edge-coupled half-wave microstrip):
    line width       = 3.058975e-3 m  (3.0590 mm)
    resonator length = 4.107003e-2 m  (41.0700 mm)
    gap[0] = 2.804932e-3 m  (2.8049 mm)   target_k = 0.069048
    gap[1] = 3.291505e-3 m  (3.2915 mm)   target_k = 0.056575
    gap[2] = 3.291505e-3 m  (3.2915 mm)   target_k = 0.056575
    gap[3] = 2.804932e-3 m  (2.8049 mm)   target_k = 0.069048
VERDICT: PASS
```

The process exits **0** on PASS and **1** on FAIL, so the command drops straight
into a CI gate or a shell `&&` chain.

Reading the report top to bottom:

- **prototype g-values** — the order-5 lowpass-prototype element values
  (`g0 … g6`), the published Chebyshev 0.5 dB N=5 values that every downstream
  stage maps from.
- **external Q** — `Qe_in` / `Qe_out`, the input/output coupling that ties the
  end resonators to the 50 Ω terminations.
- **coupling matrix M** — the normalized symmetric tridiagonal coupled-resonator
  matrix (zero diagonal ⇒ synchronously tuned).
- **mask** — the achieved in-band ripple and return loss versus the spec values,
  then one line per stopband point with the achieved rejection, the requirement,
  and `OK` / `UNDER`.
- **physical dimensions** — the edge-coupled microstrip line width, resonator
  length, and the inter-resonator gaps that realize each target coupling on the
  chosen substrate. If a coupling cannot be realized on that substrate, the
  command prints a diagnostic and exits non-zero (it never silently skips the
  dimensions).
- **VERDICT** — `PASS` iff in-band ripple ≤ spec, in-band RL ≥ spec, and every
  stopband point meets its rejection.

> **Note on scope.** The CLI `synth` subcommand always dimensions the
> **edge-coupled** realization. To explore the other five techniques
> (hairpin / combline / interdigital / lumped / stepped-impedance), use the
> studio, which routes each into its own dimensioner.

### Render the response with the spec mask

Add `--plot` to also draw `|S21|` (dB) with the mask's forbidden regions shaded —
a passband floor at `−passband_ripple_db` over the band edges, and a ceiling at
`−rejection` over a ±2 % band around each stopband point:

```bash
./target/release/yee filter synth \
    crates/yee-cli/tests/fixtures/cheb_bpf.toml \
    --output cheb_bpf.s2p --plot cheb_bpf.png
```

`--output` overrides the Touchstone path; `--plot` takes a `.png` or `.svg`
target (the extension picks the format). A passing design keeps its `|S21|`
trace out of every shaded box.

### Export a board: Gerber and KiCad

To get fabrication files for the dimensioned edge-coupled board, add `--gerber`
and/or `--kicad-pcb` (and `--layout-svg` for a preview):

```bash
./target/release/yee filter synth \
    crates/yee-cli/tests/fixtures/cheb_bpf.toml \
    --gerber cheb_bpf.gbr --kicad-pcb cheb_bpf.kicad_pcb --layout-svg cheb_bpf.svg
```

All three are emitted from the *same* layout (built once via
`dimension_edge_coupled_layout`), so the SVG preview, the Gerber, and the KiCad
board can never diverge. The KiCad `.kicad_pcb` opens directly in KiCad 7+; the
Gerber carries the copper (`F.Cu`) and board outline (`Edge.Cuts`) layers.

To run the same synthesis as a registered validation gate rather than a one-off,
`yee validate synth` runs the Filter Phase F0 gates (`synth-*` / `filt-*`).

## In the studio

Open the live studio at **<https://vsxh2001.github.io/Yee/studio/>** (or run it
locally with `dx serve` per the prerequisites). It opens seeded with a
satisfiable Chebyshev bandpass design and presents the stage rail down the left
edge. Everything re-derives live as you edit — there is no server, no build step:
the whole engine runs in your browser as WebAssembly.

Walk the stages:

1. **Spec** — a live editable form: frequency, fractional bandwidth, order,
   ripple, return loss, stopband points; pick Chebyshev or Butterworth. Every
   edit re-runs synthesis and re-grades against the mask.
2. **Technique** — the recommendation (primary + rationale + ranked
   alternatives) and the topology gallery. **Compare** the techniques' graded
   metrics in a table and overlay their `|S21|` traces on the spec mask, then
   pick one. Selecting a technique routes the rest of the rail into that flow.
3. **Synthesis** — the prototype g-values, the coupling matrix grid (or LC
   ladder for the lumped flow), the external Q, and the swept ideal
   `|S21|` / `|S11|` vs the spec mask, with a coloured **PASS** / **FAIL**
   verdict.
4. **Components** + **Tolerance** *(lumped flow only)* — the E24/E96 bill of
   materials (`select_components`) and a Monte-Carlo **yield** estimate
   (`monte_carlo_yield`) that perturbs each part by its series tolerance and
   reports the fraction of builds that still pass the mask.
5. **Layout** — the parametric microstrip (or lumped) board, rendered from the
   dimensioner.
6. **Verify** — the active flow's **circuit-level** mask metrics (see below).
7. **Export** — a parameter sheet plus **Gerber**, **KiCad**, **Touchstone**, and
   BOM downloads, generated client-side.

## Verify is circuit-level — full-wave EM is the deferred frontier

This is the one claim to be precise about. The studio's **Verify** stage grades
your design at the **circuit / synthesis level**, not with a full-wave EM solve
of the actual board:

- For the **lumped** flow, Verify grades the **realized LC ladder** (the actual
  E24/E96 part values, `ladder_s21`) against the spec mask — a genuine
  circuit-level check that the *built* ladder meets spec
  (`VerifyLevel::RealizedLadder`).
- For the **distributed** band-pass and stepped-impedance low-pass flows, Verify
  grades the **synthesized ideal** coupled-resonator response against the mask
  (`VerifyLevel::SynthesizedIdeal`).

What Verify does **not** do is simulate the physical microstrip board with a
full-wave EM solver — metal thickness, dielectric loss, dispersion, and the real
inter-resonator coupling are *not* captured. Full-wave EM verification of a
fabricated planar filter is a **deferred research frontier** for Yee (ADR-0133):
the FDTD path to a board-level S₂₁ currently floors at a single-cell aperture
port, and the route past it (a multi-cell aperture port or an FEM driven-sweep)
is a multi-week effort that has not been built. The studio states this honestly
in the Verify panel and never fabricates EM numbers; this tutorial does the
same. The Touchstone `.s2p` you export is the **ideal closed-form** response, not
an EM result.

## What you get, and what's next

From this flow you walk away with:

- a Touchstone **`.s2p`** of the ideal response — open it in any S-parameter
  tool (or feed it back through Yee's `yee-io` Touchstone reader);
- a spec-mask **PASS/FAIL** verdict — a machine-checkable gate on whether the
  synthesized response meets your requirements;
- physical **dimensions** + a parametric **board layout** for your chosen
  technique;
- (lumped flow) a **BOM** and a Monte-Carlo **yield** estimate;
- **Gerber** + **KiCad** fabrication files.

The open frontier is the **full-wave EM verification** of the fabricated board —
closing the ideal-vs-physical gap that Verify is honest about today. For the EM
solvers that work will lean on, see
[Tutorial 3 — FDTD cavity resonance](03-fdtd-cavity.md) and
[Tutorial 16 — FDTD lumped LC resonance from
Python](16-fdtd-lumped-lc-resonance-from-python.md). For the natural-language
route into a spec, see [Tutorial — Natural-language design
surface](04-nl-design-surface.md).
