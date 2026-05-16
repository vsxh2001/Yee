# patch-2g4

A 2.4 GHz rectangular patch antenna on FR-4 — the Phase 0 → 1 demo for
`yee-mom`. The example hand-builds a small planar mesh (~50 triangles) of
a 29.2 mm × 38.0 mm copper patch (tag `1`) over a ground apron (tag `2`),
with the port edge identified by the shared boundary between the two
tagged regions. It then runs a 21-point frequency sweep from 2.0 GHz to
3.0 GHz through `PlanarMoM::default().run(...)`.

The target geometry resonates at ~2.4 GHz when fed against a standard
FR-4 substrate (εr ≈ 4.4, tan δ ≈ 0.02, h ≈ 1.6 mm). However, **the
substrate is not in the mesh** — `yee-mom` models it through a multilayer
Green's function that lands in Phase 1.1. Until then, `PlanarMoM::run`
returns `Error::Unimplemented` and the example exits 0 with a clear
message pointing at the dependency.

## Run

```bash
cargo run --release --bin patch-2g4
```

## Dependencies on upcoming work

- **Phase 1.0** (Track A — `feature/phase-1-0-mom-dipole`): basic
  free-space MoM kernel. Once merged the solver stops returning
  `Unimplemented` for the simple-mesh path.
- **Phase 1.1**: multilayer dielectric Green's functions (DCIM / SDP).
  Required for the patch's resonance to land at the expected ~2.4 GHz on
  the FR-4 stack-up. Until then the free-space-only solver will be
  *qualitatively* wrong for this geometry.

## Expected output today (Phase 0)

```
patch-2g4: building rectangular patch mesh (29.2 mm × 38.0 mm, NX=4, NY=5)
patch-2g4: mesh has 34 vertices, 48 triangles (40 patch / 8 ground)
patch-2g4: invoking PlanarMoM::run over [2.00, 3.00] GHz, 21 points
patch-2g4: PlanarMoM::run is a Phase 0 stub (PlanarMoM::run not implemented in phase 0).
patch-2g4: accurate results require Phase 1.1 multilayer Green's functions (FR-4 εr ≈ 4.4, h ≈ 1.6 mm). Re-run once that lands.
patch-2g4: done.
```

## Expected output once Phase 1.1 lands

A Touchstone `.s1p` written to `target/example-output/patch-2g4.s1p`
containing `S11(f)` over 2.0–3.0 GHz, with a clear resonance near 2.4 GHz
(`|S11|` minimum, return loss > 10 dB at resonance).
