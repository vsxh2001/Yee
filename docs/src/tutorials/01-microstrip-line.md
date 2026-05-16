# Tutorial 1 — Hello, microstrip

This is the shortest end-to-end walk through Yee that exists. You will
clone the repo, build the workspace, run the bundled `microstrip-line`
example, and read the source it just executed. By the end you should
understand the *shape* of the meshing pipeline — even though the actual
Gmsh-backed mesher is gated behind a feature flag that Phase 1.mesh.0
will turn on. That gating is deliberate and is part of what this tutorial
explains.

## Goal

Run `cargo run --release --bin microstrip-line` from a freshly cloned
workspace, see a clean exit, and read the source to learn what the
binary *would* do once the `gmsh` feature is wired in. The example
exists today as a typed compile-time smoke test of the pipeline
`Session::new -> import_step -> mesh -> tris`. We are not yet generating
real triangles. We are confirming that the API surface compiles, links,
and runs.

## Prerequisites

- **Rust 1.88 or newer.** Older toolchains will fail to compile the
  workspace; the project pins its MSRV explicitly. If `rustup show`
  reports an older default, run `rustup update stable`.
- A C/C++ toolchain (most distributions install this with
  `build-essential` / Xcode CLT). Required transitively by some
  workspace deps.
- **Optional: Gmsh SDK 4.15 or newer.** Set `GMSH_SDK_ROOT` to the SDK
  install prefix to exercise the real mesh path once Phase 1.mesh.0
  ships. You can complete this tutorial without it.

You do not need Python, CUDA, or any plotting backend for this tutorial.

## Clone and build

```bash
git clone https://github.com/yee-em/yee.git
cd yee
cargo build --release
```

The first build pulls a fair number of dependencies (nalgebra, ndarray,
clap, etc.) and takes a few minutes. Subsequent builds are incremental.

## Run

```bash
cargo run --release --bin microstrip-line
```

Expected output today (Phase 0):

```
microstrip-line: attempting to open a yee-mesh Session
microstrip-line: yee-mesh built without the `gmsh` feature; mesh pipeline is gated until Phase 1.mesh.0.
microstrip-line: rebuild with `--features yee-mesh/gmsh` (Phase 1.mesh.0) to exercise the real pipeline.
microstrip-line: done.
```

The process exits with code 0. That is the intended outcome. The binary
asked `yee-mesh` to construct a `Session`, the crate replied
`Error::NotEnabled` because the `gmsh` cargo feature is off in the
default build, and the example printed a clear explanation rather than
panicking.

If you supply a STEP file path, the binary acknowledges it but still
exits cleanly:

```bash
cargo run --release --bin microstrip-line -- --step path/to/board.step
```

## What the example actually does

Open `examples/microstrip-line/src/main.rs`. The entire driver is about
eighty lines and is worth reading top to bottom. The interesting calls
are:

```rust
match yee_mesh::Session::new() {
    Ok(mut session) => {
        session.import_step(&step)?;
        session.mesh(2)?;
        let tris = session.tris()?;
        println!("microstrip-line: mesh has {} vertices, {} triangles",
                 tris.vertices.len(), tris.n_tris());
    }
    Err(yee_mesh::Error::NotEnabled) => { /* explain and exit 0 */ }
    Err(other) => return Err(anyhow::anyhow!("unexpected: {other}")),
}
```

This is the canonical Yee meshing pipeline: open a session, import a
geometry file (STEP today; IGES and KiCad PCBs later), call `mesh(dim)`
to triangulate the surface (`dim = 2`), and pull the triangle list out
as a `TriMesh`. The same `TriMesh` is what `yee-mom` consumes downstream
(see the next tutorial), so this pipeline is the upstream half of every
planar-MoM simulation in the studio.

The three-armed `match` is what makes the example robust:

1. **`Ok(session)`** — the `gmsh` feature is on and meshing succeeds.
   The example prints triangle counts. (Not reachable in today's
   default build.)
2. **`Err(Error::NotEnabled)`** — the `gmsh` feature is off. The
   example prints a friendly message and exits 0. *This is what you
   saw above.*
3. **Any other error** — propagated through `anyhow` and surfaced.

## Enabling Gmsh (optional, advanced)

Once Phase 1.mesh.0 lands the bindgen-generated FFI against `gmshc.h`,
the build is:

```bash
export GMSH_SDK_ROOT=/path/to/gmsh-4.15-sdk
cargo run --release --bin microstrip-line \
    --features yee-mesh/gmsh -- --step examples/microstrip-line/reference/microstrip.step
```

You will then see a line like `mesh has N vertices, M triangles` and
the same binary becomes a useful sanity check for any STEP file you
feed it. Today, that build path is unreachable — `Session::new` is a
`todo!()` placeholder behind the feature flag. The README in
`examples/microstrip-line/README.md` tracks the current status.

## Why this matters

Most simulation studios hide the mesh stage behind a GUI or a giant
configuration file. Yee exposes it as a small typed API surface that
you can compile against today and that will silently start doing real
work once the FFI ships. The tradeoff is that early adopters see a
"feature not enabled" message in places where a finished tool would
print triangle counts; the upside is that nothing about your driver
code has to change when Phase 1.mesh.0 merges.

## Next

Move on to [Tutorial 2 — Half-wave dipole from
Python](02-dipole-from-python.md), where the solver side of the pipeline
*is* implemented and you can see real impedance numbers come back from a
hand-built cylinder mesh.
