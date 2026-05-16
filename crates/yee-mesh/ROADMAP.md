# yee-mesh — Roadmap

## Phase 0 (months 0–6)
- [ ] `TriMesh` struct + tag-vector + basic invariants
- [ ] `build.rs` invokes `bindgen` against `gmshc.h` from $GMSH_SDK
- [ ] Safe wrappers: `Session::new` / `drop`, `import_step`, `mesh(dim)`, `tris()`
- [ ] Example: `examples/microstrip-line.rs` — extrude trace + GND, mesh, count tris
- [ ] CLI smoke test invoked from `yee-cli validate mesh`

## Phase 1 (months 6–18)
- [ ] Size fields (constant, distance-from-curve, anisotropic)
- [ ] Physical group tagging API
- [ ] KiCad → OCC stack-up importer
- [ ] Via / slot handling
- [ ] Adaptive remesh hooks for solver-driven refinement

## Phase 2 (months 18–30)
- [ ] Volumetric Cartesian / hex / tet meshes
- [ ] Non-uniform Cartesian grids with stability validation
- [ ] Voxelizer (OCC shape → Yee grid occupancy)
- [ ] Conformal Dey-Mittra fractional-cell support

## Phase 4
- [ ] Edge-element conformal meshes for FEM

## Validation gates per phase
- Phase 0: importing a simple STEP cube yields the analytical surface-triangle count within ±5% across refinement levels.
- Phase 1: KiCad import of a published test PCB matches Gerber outline within ±10 µm.
- Phase 2: voxelization conserves volume to ±1 cell-volume on a unit sphere.

## Watch-outs
- Gmsh SDK install path must be discoverable: `GMSH_SDK_ROOT` env var supported in `build.rs`.
- C API surface is large (~500 functions); we expose only what we need and grow the surface deliberately.
