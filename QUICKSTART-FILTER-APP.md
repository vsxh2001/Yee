# Yee Filter-Design — Quickstart

The end-to-end RF filter-design flow that ships today: a spec goes in, and a
synthesized response, physical microstrip dimensions, a layout preview, and
manufacturing files (Gerber + KiCad) come out — from the CLI or the in-browser
Studio. (Full-wave FDTD verification of the realized coupling is landing as
F1.1b.1, ADR-0108.)

This is a walking-skeleton product: the **planar edge-coupled** topology on a
single substrate is the path wired end to end. See `FILTER-DESIGN-ROADMAP.md`
for the broader plan (waveguide / lumped tracks, more topologies, surrogate-BO
EM-in-the-loop).

---

## 1. CLI: spec → response + dimensions + manufacturing files

```bash
# Synthesize, sweep the ideal response, dimension onto FR-4, and emit every
# downstream artifact in one shot. All --flags are optional and share ONE
# computed layout (the SVG, Gerber, and KiCad board can never diverge).
yee filter synth my_filter.toml \
  --output    my_filter.s2p          \  # Touchstone S-parameters (ideal response)
  --plot      my_filter.png          \  # |S21| with the spec-mask overlaid
  --eps-r     4.4                    \  # substrate epsilon_r   (default FR-4 4.4)
  --h-mm      1.6                    \  # substrate height (mm)  (default FR-4 1.6)
  --layout-svg my_filter.svg         \  # top-view geometry preview
  --gerber    my_filter.gbr          \  # single-copper-layer RS-274X Gerber
  --kicad-pcb my_filter.kicad_pcb       # KiCad 7 board (opens in the PCB editor)
```

Exit code is `0` if the synthesized response PASSES the spec mask, `1` on FAIL,
so `yee filter synth` doubles as a CI-able design check. The printed block shows
the realized line width, resonator length, and inter-resonator gaps.

A minimal `my_filter.toml` (Chebyshev band-pass) — see
`docs/src/tutorials/` for the full tutorial:

```toml
response       = "Bandpass"
f0_hz          = 2.0e9
fbw            = 0.10
order          = 5
z0_ohm         = 50.0
[approximation]
Chebyshev = { ripple_db = 0.5 }
[mask]
passband_ripple_db = 0.5
return_loss_db     = 10.0
stopband           = [[2.4e9, 30.0]]
```

### Board outline + opening in KiCad

The `.kicad_pcb` carries the copper traces on `F.Cu` and a board outline on
`Edge.Cuts`; open it directly in KiCad's PCB editor. For a fab hand-off, the
Gerber copper (`--gerber`) plus the Edge.Cuts board outline
(`yee_export::layout_to_gerber_outline`) are the RS-274X equivalents.

---

## 2. Validation dashboard

```bash
yee validate --list      # every gate, incl. the filter pipeline:
                         #   synth-001/002, filt-001  (synthesis)
                         #   coupled-001              (coupled-line model vs Steer)
                         #   dim-001                  (dimensional-synthesis round-trip)
                         #   gerber-001               (Gerber RS-274X structure)
```

These filter gates are pure-math/text (millisecond-scale) and run without any
EM solve.

---

## 3. Web UI (in-browser Studio)

The same `yee-studio` egui app runs natively and in the browser (one codebase,
ADR-0089). To run the web build locally:

```bash
rustup target add wasm32-unknown-unknown
cargo install --locked trunk
trunk serve crates/yee-studio/index.html
# → open the printed http://127.0.0.1:8080
```

The Studio gives live synthesis, editable substrate (ε_r / h), the dimensioned
geometry, and a top-view layout canvas. CI builds the deployable static bundle
(`crates/yee-studio/dist`) on every push (the `wasm-build` job uploads it as the
`yee-studio-web` artifact).

The native desktop app is just `cargo run -p yee-studio`.

---

## 4. What runs the EM (status)

The flow above is closed-form (synthesis + analytic coupled-line models) — fast,
no solver. The first **full-wave** step, `yee-voxel::run_coupled_pair` (FDTD
coupled-resonator coupling extraction, F1.1b.1 / ADR-0108), validates that the
dimensioned geometry realizes the target coupling; its multi-minute FDTD gate
runs in CI (`fdtd-coupling-gate`). Surrogate-BO EM-in-the-loop refinement
(F1.2.1) builds on it.
