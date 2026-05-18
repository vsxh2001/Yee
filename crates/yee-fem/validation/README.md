# yee-fem — Validation

Every solver feature is held against a canonical published benchmark
before it ships. No exceptions (CLAUDE.md §4). This page tracks the
Phase 4 FEM eigenmode validation rollup.

## Canonical references

- Pozar, *Microwave Engineering* (4th ed.) — §6.3 rectangular metallic
  cavity, analytic TE_{mnp} / TM_{mnp} frequencies (eq. 6.42), wall-loss Q.
- Jin, *The Finite Element Method in Electromagnetics* (3rd ed.) —
  Ch. 9 (Nedelec edge elements on tetrahedra), Ch. 10 (eigenvalue
  problems and cavity resonators).
- Petosa, *Dielectric Resonator Antenna Handbook* — cylindrical DRA
  validation case for the deferred `fem-eig-003` (Phase 4.fem.eig.3).

## Cases — Phase 4.fem.eig.0 (walking skeleton)

| ID | Case | Tolerance | Wall-time |
|----|------|-----------|-----------|
| `fem-eig-001 (TE_{101})` | WR-90-based rectangular metallic cavity, `a = 22.86 mm`, `b = 10.16 mm`, `d = 30 mm`, lossless air fill, `(nx, ny, nz) = (12, 9, 15)` Kuhn 6-tet brick mesh (9720 tets, 10 107 interior DoFs); lowest mode vs the Pozar §6.3 eq. 6.42 analytic TE_{101} | `|f_FEM − f_TE101_analytic| / f_TE101_analytic ≤ 0.3 %` | `< 60 s` in `--release` (informational; ~10 s observed) |
| `fem-eig-001 (mode-10 ordering)` | Same geometry / mesh / solve as the TE_{101} row; the ten lowest measured eigen-frequencies are compared mode-by-mode against the Pozar §6.3 analytic TE/TM table (eq. 6.42), with the no-spurious-mode-below-TE_{101} sanity check enforced post-solve | ±1 % pairwise per mode; every returned eigenvalue `> 0.5 · k₀_TE101²` so no gradient-cluster mode appears below TE_{101} | covered by the same test |

Both rows are exercised by `crates/yee-validation/tests/fem_eig_001_rectangular_cavity.rs`,
which drives the public `yee_validation::run_fem_eig_001_rectangular_cavity`
helper.

### Findings surfaced during the T7 landing

* **Pozar TE_{101} reference value vs cavity dimensions.** Phase 4
  spec §9 and the T7 agent brief both cite `f_{101} ≈ 9.660 GHz` as
  the TE_{101} target for `(a, b, d) = (22.86 mm, 10.16 mm, 30 mm)`.
  That value is the Pozar 4th ed. *worked example* for a WR-90 cavity
  with `d = 20 mm`, not `d = 30 mm`. Applying Pozar eq. 6.42 directly
  to the spec's stated dimensions yields `f_{101} ≈ 8.244 GHz`. The
  mode-10 ordering gate (2) also requires consistency with the Pozar
  table evaluated at the same `(a, b, d)`, so the only self-consistent
  reference is the formula. The driver computes the TE_{101} target
  inline (`fem_eig_001_f_te101_hz`) rather than hardcoding `9.660e9`.
  Surfaced as a finding rather than fixed in the spec — that lane
  belongs to a future spec-edit track.

* **Mesh resolution requirement.** The T7 brief calls for the default
  gate to run at `(nx, ny, nz) = (8, 6, 10)` (2880 tets). At that
  resolution the TE_{101} bound (gate 1) is met comfortably (0.19 %
  measured error) but several mode-10 modes whose field profile varies
  across the narrow `b = 10.16 mm` direction land 1.2 %–1.4 % low,
  exceeding the ±1 % bound on gate (2). Per the brief's escape hatch
  ("refine to (12, 9, 15) and retry"), we run the default gate at
  `(12, 9, 15)` (9720 tets), where every mode-10 mode lands within
  ±0.6 %. Wall-time at the refined resolution is ~5–10 s in
  `--release`, well inside the 60 s informational and 5 min
  `#[ignore]` thresholds — no fallback is needed.

* **Shift `σ` choice for the Phase 4 T5 escape-hatch eigensolver.**
  The brief specifies `σ = 0.5 · k₀_TE101²` ("below the smallest
  physical mode but above the gradient-kernel cluster at 0"). That
  literal value is an unfortunate boundary case for the deflated
  inverse-power iteration that ships as the T5 escape-hatch impl
  (`yee_fem::InverseIterEigen`): TE_{101} at `k² ≈ 2σ` and the
  gradient kernel at `k² = 0` produce identical `|θ|` magnitudes, so
  inverse-iteration has no preference for the physical mode. We lift
  `σ` to `2.5 · k₀_TE101²` (sits between the 8th and 9th physical
  modes of the Pozar table) so the gradient cluster is decisively
  outranked and every requested mode converges in ascending `k²` order.
  This dependency on the shift heuristic is a documented limitation
  of the escape-hatch impl; the spec §8 `SparseEigen` trait abstracts
  the solve, so a future LOBPCG / ARPACK swap removes the heuristic
  in one PR.

## Deferred cases

- `fem-eig-002` (lossy-cavity Q-factor): Phase 4.fem.eig.2. Requires
  complex `ε_r` end-to-end and a complex generalized eigensolve;
  validated against Pozar §6.3 wall-loss Q. Out of scope at v0.
- `fem-eig-003` (cylindrical DRA): Phase 4.fem.eig.3. Petosa DRA
  Handbook ch. 3 tabulation. Out of scope at v0.

## Running

```bash
# Default-features production gate.
cargo test -p yee-validation --release --test fem_eig_001_rectangular_cavity
```

Driver returns a `yee_validation::FemEigValidationResult` carrying the
measured + analytic frequencies, per-mode relative errors, headline
TE_{101} bound, mode-10 RMS error, and pass/fail status.

## Mesh-quality note

The Kuhn 6-tet decomposition of an axis-aligned brick (`yee-mesh`'s
`TetMesh3D::cavity_uniform`) is well-conditioned by construction: each
brick yields six congruent orthoschemes with identical signed volume,
so the lowest mode hits the Pozar 4-significant-digit table comfortably
inside the ±0.3 % bound on the default `(8, 6, 10)` mesh. Arbitrary
Gmsh-imported tet meshes are deferred to Phase 4.fem.eig.1 along with a
`MeshSizeFactor ≤ 0.05·λ` recommendation; v0 uses the hand-rolled
cavity only (Phase 4 spec §11 risk register).

## Plot artifacts

Phase 4.fem.eig.0 ships without plot artifacts — the validation gate
exercises the eigen-frequency table only. Mode-profile visualisation is
the optional T9 deliverable (mdBook tutorial `04-fem-cavity-eigenmode.md`)
that lands after the production gate.
