# yee-io — Roadmap

## Phase 0 (months 0–6)
- [ ] Touchstone v1.1 grammar (option line, comments, frequency units MHz/GHz/Hz)
- [ ] Reader: produces `touchstone::File` for `.s1p` through `.s4p`
- [ ] Writer: stable line ordering, deterministic output
- [ ] Round-trip property test: read → write → read = identical struct
- [ ] Corpus of published reference files in `validation/fixtures/touchstone/`

## Phase 1 (months 6–18)
- [ ] Generic `.sNp` (n > 4) support
- [ ] Touchstone v2 syntax (option list, frequency-dependent Z₀, mixed-mode)
- [ ] OCC STEP / IGES import via `opencascade-rs` 0.2+
- [ ] STL import for free-form bodies
- [ ] KiCad PCB import (stack-up, layers, vias, drills)
- [ ] HDF5 output via `hdf5` crate for field arrays
- [ ] Arrow IPC writer for surrogate training datasets

## Phase 2 (months 18–30)
- [ ] Gerber import for legacy layers
- [ ] DXF / SVG outlines for antenna teaching examples
- [ ] Rerun-sdk helpers (stream solver internals to the Rerun viewer)

## Validation gates per phase
- Phase 0: 100% round-trip on the curated reference corpus; lint-clean S-parameter passivity check on read.
- Phase 1: OCC import of every test geometry in `validation/cad/` matches OCC's own surface-area to ±0.1%.
- Phase 1: KiCad import of a published demo PCB matches Gerber outline to ±10 µm.

## Risks
- Touchstone has many vendor dialects (HFSS adds non-standard comments; Sonnet has its own header conventions). We parse permissively and emit strictly.
- OCCT 8.0 lands ~May 2026; expect `opencascade-rs` API churn after that. Pin and migrate deliberately.
