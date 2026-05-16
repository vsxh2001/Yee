# half-wave-dipole

Phase 0 walking-skeleton example for the planar Method of Moments solver
(`yee-mom`). The example constructs a minimal two-triangle planar mesh that
represents a free-space half-wave dipole at 300 MHz, with the port edge
identified by the convention "shared edge between two differently-tagged
triangles." It then hands the mesh to `PlanarMoM::default().run(&mesh, freq)`
and writes a Touchstone `.s1p` containing the resulting S-parameters.

Until Track A (`feature/phase-1-0-mom-dipole`) lands real physics, the
solver is still a Phase 0 stub: it returns `Error::Unimplemented`. This
example handles that gracefully — it prints the unimplemented message and
writes a clearly-labelled placeholder `.s1p` to `target/example-output/`,
exiting `0`. Once the MoM kernel ships, the same driver will compute the
classical free-space dipole impedance `Z ≈ 73 + j42 Ω`.

## Run

```bash
cargo run --release --bin half-wave-dipole
```

Output is written to `target/example-output/dipole.s1p` (gitignored).

## Expected output today (Phase 0)

```
half-wave-dipole: building two-triangle planar mesh
half-wave-dipole: mesh has 5 vertices, 2 triangles, tags = [1, 2]
half-wave-dipole: invoking PlanarMoM::run at 300.000 MHz
half-wave-dipole: PlanarMoM::run is a Phase 0 stub (PlanarMoM::run not implemented in phase 0).
half-wave-dipole: writing placeholder S-parameters until Track A merges.
half-wave-dipole: wrote target/example-output/dipole.s1p
half-wave-dipole: done.
```

## Expected output once Phase 1.0 lands

A `.s1p` file containing `S11(300 MHz)` for a free-space dipole consistent
with the closed-form `Z ≈ 73 + j42 Ω` (matched to `Z0 = 50 Ω`, so
`|S11| ≈ 0.38` with a small positive phase).
