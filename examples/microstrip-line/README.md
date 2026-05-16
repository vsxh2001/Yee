# microstrip-line

Phase 0 smoke example for the meshing pipeline (`yee-mesh`). The binary
attempts to open a `yee_mesh::Session`, optionally import a STEP file, and
report the triangle count of the resulting surface mesh.

Today (Phase 0), `yee-mesh` ships without the `gmsh` feature by default —
the bindgen FFI to the Gmsh SDK is gated behind that feature and is the
deliverable of **Phase 1.mesh.0**. In the default build, `Session::new()`
returns `Error::NotEnabled` and this example prints a clear message before
exiting 0. Once the SDK is wired in, the same driver will exercise the
import → mesh → tris pipeline end-to-end.

## Run

```bash
# Default: prints a "gmsh feature not enabled" message and exits 0.
cargo run --release --bin microstrip-line

# With a STEP file (only meaningful once Phase 1.mesh.0 ships):
cargo run --release --bin microstrip-line -- --step path/to/microstrip.step
```

## Dependencies on upcoming work

- **Phase 1.mesh.0** ships the bindgen-generated FFI against `gmshc.h`
  and unblocks the `--features yee-mesh/gmsh` build. Until then the
  example serves as a typed compile-time smoke test of the surface
  `Session::new → import_step → mesh → tris`.

## Expected output today (Phase 0)

```
microstrip-line: attempting to open a yee-mesh Session
microstrip-line: yee-mesh built without the `gmsh` feature; mesh pipeline is gated until Phase 1.mesh.0.
microstrip-line: rebuild with `--features yee-mesh/gmsh` (Phase 1.mesh.0) to exercise the real pipeline.
microstrip-line: done.
```
