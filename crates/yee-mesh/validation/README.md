# yee-mesh — Validation

Gmsh-gated. CI runs the data-structure tests on every PR; the FFI tests require a Gmsh SDK and run on the meshing runner.

## Cases

| ID | Description | Reference | Tolerance | Phase |
|----|-------------|-----------|-----------|-------|
| `mesh-001` | TriMesh round-trip: insert tris, read back, count | analytical | bit-exact | 0 |
| `mesh-002` | Gmsh import unit cube STEP → ~12 tris @ coarse | OCC reference | exact at coarsest | 0 |
| `mesh-003` | Microstrip-line example builds + meshes without error | exit 0 | — | 0 |
| `mesh-004` | KiCad PCB → outline boundary recovery | Gerber edge | ±10 µm | 1 |
| `mesh-005` | Voxelization of unit sphere conserves volume | 4π/3 | ±1 voxel | 2 |
| `mesh-006` | Refinement convergence — uniform refinement halves edge length | doubling | ±2% | 1 |

## Running

```bash
# Data-only tests (no Gmsh required)
cargo test -p yee-mesh

# Full FFI tests (Gmsh SDK on host)
GMSH_SDK_ROOT=/opt/gmsh-sdk cargo test -p yee-mesh --features gmsh
```

## CI

- PR: data-only build + tests.
- Nightly: full FFI on a runner with the Gmsh 4.15+ SDK preinstalled.
