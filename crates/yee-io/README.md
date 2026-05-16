# yee-io

> Format I/O for Yee: Touchstone, CAD (STEP/IGES/KiCad), HDF5/Arrow.

## Scope

### Phase 0
- **Touchstone v1.1** reader/writer: `.s1p`, `.s2p`, `.s3p`, `.s4p` minimum; generic `.sNp` support
- Round-trip preservation of reference impedance, frequency unit, parameter type (S/Y/Z), format (RI/MA/DB)

### Phase 1
- **OpenCascade CAD ingest** via `opencascade-rs` 0.2+: STEP, IGES, STL, DXF, SVG
- **KiCad PCB import**: stack-up, copper layers, vias, drills
- **Touchstone v2** (frequency-dependent reference impedance, mixed-mode S-parameters)
- **HDF5 output** for field arrays (E, H, J) at sweep points
- **Arrow IPC** for surrogate training datasets

### Phase 2
- **Gerber import** for arbitrary copper layers
- **`rerun` SDK** sink helper crate (already optional via workspace dep)

## Feature flags

| Flag | Effect |
|------|--------|
| `touchstone` | default. Touchstone v1.1 reader/writer. |
| `opencascade` | Pull in `opencascade-rs` for CAD. 5–15 min cold compile (OCCT submodule). |

CI caches OCCT build artifacts aggressively (`sccache`).

## Validation

See [`validation/README.md`](validation/README.md). Touchstone round-trips against a corpus of published reference files (Sonnet examples, ADS sample files, vendor reference data).

## Roadmap

See [`ROADMAP.md`](ROADMAP.md).
