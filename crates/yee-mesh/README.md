# yee-mesh

> Gmsh-backed meshing layer for Yee.

## Why in-tree FFI

`rgmsh` / `gmsh-sys` (mxxo) — last updated November 2019, targets Gmsh 4.4.1, six years stale. **Unmaintained.** We regenerate bindings ourselves with `bindgen` against the current `gmshc.h` and write a thin safe wrapper here. That gives us:

- Control over the Gmsh version we link against (4.15+ as of 2026, watching 4.16).
- A safe API that maps cleanly onto `yee-core`'s mesh trait (Phase 1).
- An escape hatch: shell out to the Gmsh CLI or its Python API when the binding lags.

## Scope

### Phase 0
- `TriMesh` data structure + tagging
- `build.rs` generating `bindings.rs` from `gmshc.h` when `--features gmsh` is on
- Safe wrapper for: initialize, finalize, OCC geometry import, mesh generation, element retrieval
- Microstrip-line example: extrude a rectangular trace + ground, mesh, dump `TriMesh`

### Phase 1
- Refinement controls (size fields, anisotropy along ports)
- Physical groups → port/material tag mapping
- KiCad PCB → Gmsh OCC stack-up importer
- Robust handling of overlapping copper, vias, slot apertures

### Phase 2
- Volumetric hex/tet meshes for FDTD voxelization
- Non-uniform Cartesian grid generator with stability fixes
- Conformal Dey-Mittra cell classification

## Feature flags

| Flag | Effect |
|------|--------|
| (none) | Data structures only. CI safe on hosts without Gmsh. |
| `gmsh` | Build + link against an installed Gmsh SDK 4.15+. |

## License notes

Gmsh is **GPL v2+ with linking exception** — FAQ-confirmed compatible with our GPL v3 host. Document the dependency in `THIRD_PARTY_LICENSES.md` at the workspace root (Phase 0 deliverable).

## Validation

See [`validation/README.md`](validation/README.md).

## Roadmap

See [`ROADMAP.md`](ROADMAP.md).
