# FDTD Subgridding — Theory of Operation

This page is the theory-of-operation reference for the FDTD
subgridding scheme planned for Phase 2.fdtd.7 of `yee-fdtd`. It is
the structural analog of [CPML / NTFF / TF-SF](./fdtd-details.md) —
interface construction, update interleaving, stability analysis,
validation — for the case where a uniform Yee grid is no longer
adequate and a locally refined region must be nested inside the
coarse parent. Audience: an engineer implementing the walking
skeleton (`SubgridRegion`, `SubgriddedSolver`), debugging a
late-time interface instability, or extending to higher refinement
ratios. Equations are plain-text ASCII — the mdBook build has no
math preprocessor (see `docs/book.toml`).

## 1. Introduction

The base FDTD chapter ([`fdtd.md`](./fdtd.md)) constructs Yee's
staggered grid with a single global cell size `dx = dy = dz`. That
works for closed cavities and free-space dipoles. It breaks down
the moment a single sub-wavelength feature forces refinement of the
entire volume: a 2.4 GHz patch on 0.508 mm RO4003C
(`λ_eff / 80` thick) demands `dx ≤ 50 µm` globally, putting a
100×100×30 mm domain at `≈ 2.4 · 10⁹` cells; slot widths of
`λ / 200` waste the same effort on surrounding free space;
conductor-corner singularities benefit super-linearly from local
refinement, and staircasing corrupts radiated patterns at the
few-percent level NTFF validation cares about.

**Subgridding** refines only where the geometry demands it: a
coarse parent grid covers the bulk of the domain, and one or more
nested **fine sub-grids** cover the thin features. The price is the
coarse / fine **interface**, where field components from two grids
that disagree on spacing must be reconciled at each leapfrog step.
The classical failure mode is **late-time instability**: a
sub-percent reciprocity error per step compounds over `O(10⁴)`
steps into exponential growth of an unphysical interface eigenmode.
Most of this chapter is about why the interface contract must be
precisely reciprocal and how the Chevalier–Okoniewski scheme
arranges that. Phase 2.fdtd.7.0 ships the **walking skeleton**: one
axis-aligned cuboidal fine region nested at fixed `r = 2` inside
one parent `YeeGrid`, with linear spatial / temporal interpolation
coarse → fine and area-averaged fine → coarse closure. Sub-phases
7.1–7.5 lift the single-nest, single-ratio, non-dispersive,
CPML/TF-SF-interior restrictions one at a time. Companion spec:
`docs/superpowers/specs/2026-05-18-phase-2-fdtd-7-subgridding-design.md`.

## 2. The CFL constraint at the heart of the problem

The base chapter §4 derives the 3-D Courant–Friedrichs–Lewy bound on
the leapfrog time step; for isotropic `dx = dy = dz = h` this is

```text
c · dt  ≤  h / sqrt(3).
```

The Courant limit is **proportional to cell size**: halving `h`
halves the maximum stable `dt`. Refining one cubic region by
`r = 2`: cell count grows by `r³ = 8`; stable `dt` shrinks by `r`,
so step count grows by `2`; per-step work in the refined region
grows by `r³ · 2 = 16`. Applying the same refinement uniformly to
the **entire** domain therefore costs `16×`; restricting it to a
fraction `f` of the volume costs `16 · f` — a near-`1/f` saving any
time `f ≪ 1`. The flip side: the time-step gap between coarse and
fine is itself the interface problem. The standard fix is **time
sub-cycling**: the fine grid steps `r` times per coarse step at
`dt_fine = dt_coarse / r`. For `r = 2` that is two fine half-steps
per coarse step; §6 describes the interleaving.

## 3. Region decomposition

Phase 2.fdtd.7.0 places one **axis-aligned cuboidal fine region**
inside the parent coarse `YeeGrid`. The fine region occupies coarse
cells `lo..hi` and is itself a `YeeGrid` sized
`(r·(hi.0-lo.0), r·(hi.1-lo.1), r·(hi.2-lo.2))` fine cells at spacing
`dx_fine = dx_coarse / r`. With `r = 2` the fine origin sits on a
coarse `E`-node and every other fine `E`-node coincides with a coarse
`E`-node along each axis. On the `+x` face the Yee staggering places
**tangential `E`** (`E_y`, `E_z`) on edges in the face plane (every
coarse `E_t` edge coincides with one fine `E_t` edge; between every
pair sits one extra fine edge with no coarse counterpart) and
**normal `H`** (`H_x`) on the face itself (one coarse `H_x` cell
covering four fine `H_x` cells in a `2 × 2` pattern). Tangential `H`
and normal `E` sit off the face plane; they are interior to either
grid and do not cross the interface. The interface coupling is
therefore between coarse and fine **tangential `E`** (coarse → fine
driver) and between coarse and fine **normal `H`** (fine → coarse
closure). Symmetry permutes the same statement onto the other five
faces.

## 4. Spatial interpolation at the interface

At each fine sub-step the outer tangential `E` of the fine grid
needs a Dirichlet value sampled from the coarse grid. Two
sub-populations of fine `E_t` edges exist on the `+x` face:
**coincident** with a coarse edge (direct copy,
`E_t_fine[2j, 2k] = E_t_coarse[j, k]`), and **between** two coarse
edges (linear interpolation). For the fine edge at `2j+1` in `y`,
sandwiched between coarse `j` and `j+1`,

```text
E_t_fine[2j+1, 2k]  =  0.5 · ( E_t_coarse[j, k]  +  E_t_coarse[j+1, k] ).
```

The diagonal fine edge `(2j+1, 2k+1)` weights its four flanking
coarse values 0.25 each — a bilinear average.

**Why linear and not cubic?** Cubic Lagrange interpolation cuts the
per-step coupling error from `O(dx_coarse²)` to `O(dx_coarse⁴)`, but
at `r = 2` the per-step error is already dominated by second-order
Yee dispersion. Cubic also broadens the stencil and complicates the
reciprocity argument in §5; both costs buy nothing at v0 scope.
Linear interpolation is the Chevalier 1997 choice; Phase 2.fdtd.7.4
revisits higher-order when `r ≥ 3` makes dispersion no longer
dominant. The spatial-interpolation operator is a sparse rectangular
matrix `R: ℝ^{N_coarse} → ℝ^{N_fine}` with one or two non-zeros per
row, each `+1`, `+1/2`, or `+1/4`. `R` is **load-bearing in §5**
where its transpose appears as the closure of the energy balance.

## 5. Fine → coarse averaging

After the fine grid has completed both sub-steps for a coarse step,
the coarse grid's field on the interface must be reconciled with the
fine grid's view of the same plane. The Chevalier 1997 scheme uses
**area-weighted averaging** of the four fine cells covering each
coarse face cell:

```text
H_n_coarse[j, k]  =  0.25 · ( H_n_fine[2j,   2k  ]
                            + H_n_fine[2j+1, 2k  ]
                            + H_n_fine[2j,   2k+1]
                            + H_n_fine[2j+1, 2k+1] ).
```

The same averaging applies symmetrically to coarse `E_t` on the
interface, overwritten by an edge-average of the two coincident fine
`E_t` edges per coarse edge.

**Reciprocity argument.** Let `R` be the spatial-interpolation
operator from §4 and `T` the averaging operator above. The discrete
energy balance across the interface closes when `T = R^T` (up to a
constant geometric scale): the average pushing fine field back onto
the coarse grid must be the **transpose** of the interpolation that
pushed coarse onto fine. Linear interpolation with one-half weights
pairs naturally with an area-average using one-quarter weights
because `(1/2)² = 1/4`. A mismatch — e.g. linear interpolation
coarse → fine paired with **nearest-neighbour** fine → coarse — is
the canonical Berenger 2003 instability source: per-step it leaks
`(R^T − T)` energy onto an unphysical interface mode, and over
`O(10⁴)` steps that grows exponentially. Chevalier 1997 §IV gives
the derivation; Berenger 2003 (T-AP 51.10) shows both failure-mode
directions. The §7 stability gate experimentally checks the
transposition; the analytical guarantee is necessary but not
sufficient against indexing bugs.

## 6. Temporal sub-cycling

At `r = 2` the fine grid steps **twice** per coarse step at
`dt_fine = dt_coarse / 2`. The interleaving has to be chosen so the
fine grid always reads a temporally consistent coarse `E_t` on its
boundary and the coarse grid reads a fine-averaged `H_n` that
incorporates both fine sub-steps. The Okoniewski 1997 pattern gives
the canonical seven-stage sequence for one coarse step `n → n + 1`:

```text
1. coarse: update_h          (advances coarse H to n + 1/2)
2. coarse: apply_cpml_h, source_h
3. fine sub-step k = 1:
     a. interpolate coarse E_t at t = n + 1/4 onto fine boundary E_t
        (linear blend: 0.75 · E_t_coarse^{n} + 0.25 · E_t_coarse^{n+1})
     b. fine: update_h        (advances fine H to n + 1/4)
     c. fine: update_e        (advances fine E to n + 1/2)
4. fine sub-step k = 2:
     a. interpolate coarse E_t at t = n + 3/4 onto fine boundary E_t
        (linear blend: 0.25 · E_t_coarse^{n} + 0.75 · E_t_coarse^{n+1})
     b. fine: update_h        (advances fine H to n + 3/4)
     c. fine: update_e        (advances fine E to n + 1)
5. coarse: update_e           (advances coarse E to n + 1)
6. coarse: apply_cpml_e, source_e
7. average_fine_h_to_coarse, overwrite_coarse_e_from_fine.
```

The blend fractions `frac ∈ {0.25, 0.75}` come from the fine
sub-step midpoints relative to the coarse interval: sub-step `k = 1`
needs driver-`E_t` at `t = n + 1/4`, weight `0.25` on the later
coarse snapshot; sub-step `k = 2` needs it at `t = n + 3/4`, weight
`0.75`. Two snapshots cached in `SubgridRegion::e_t_snapshots`
supply both blends. Coarse `update_h` precedes the fine sub-steps
but coarse `update_e` follows them: stage 7's
`average_fine_h_to_coarse` closes the loop for the **next** coarse
step. The seven-stage sequence reuses the Step-1 refactor helpers
(`update_h_only`, `apply_cpml_h`, `update_e_only`, `apply_cpml_e`)
so no CPML wiring is re-implemented inside the subgrid module.

## 7. Stability and reciprocity

Late-time instability is the classical subgridding failure mode. The
mechanism is the **asymmetric-coupling eigenmode**: if coarse → fine
injects energy at a different rate than fine → coarse extracts it,
a discrete-interface eigenmode exists whose amplification factor
exceeds unity by `O(10⁻⁶)` per step — below floating-point noise per
step, compounding over `10⁴` steps to `exp(10⁻⁶ · 10⁴) ≈ 2.7×`.

The energy-conservation requirement from §5 — `T = R^T` up to
geometric scale — is the **analytical** stability condition. The
Chevalier area-average is the explicit construction of `T` for
`r = 2` linear interpolation; with that construction, the discrete
energy integral

```text
W(t)  :=  ∫ [ ε_0 |E(t)|²  +  μ_0 |H(t)|² ] dV
```

is conserved to second order in `dx_coarse` over the time-step. The
**experimental** stability gate is the round-trip energy test in
implementation-plan Step 6: initialise a Gaussian-modulated
sinusoid in the fine region, propagate forward across the interface,
reflect off PEC walls on outer faces, propagate back through the
interface, and integrate `W(t)` at `t = 0` and `t = 10⁴ · dt_coarse`
(tolerance `0.5%`). Berenger 2003 §IV is the canonical analysis;
Berenger 2006 (T-AP 54.12) is a stronger-stability Huygens-surface
fallback, escape-hatched into Phase 2.fdtd.7.x if Chevalier proves
insufficient.

## 8. CPML interaction (deferred)

The convolutional PML — see [`fdtd-details.md`](./fdtd-details.md)
§2 — adds an auxiliary state vector `Ψ` to every cell inside the PML
shell, evolving as a one-tap recursive convolution of field history
with an exponential kernel from the local stretching parameters
`(σ, κ, α)`. The CPML update is a function of field **history**, and
the history at an interface between two grids with different `dt`
and `dx` is precisely what the bare subgrid reciprocity argument is
silent on. The Chevalier–Okoniewski energy balance closes for the
bare Yee stencil but does **not** automatically close for CPML: the
discrete energy in the polarisation reservoirs leaks through the
interface in a direction the bare-grid transpose cannot cancel.
Phase 2.fdtd.7.0 therefore **forbids** co-location of a fine region
with a CPML face — `SubgridRegion::new` errors if `lo`/`hi` overlap
the parent's CPML thickness — and defers the construction to
**Phase 2.fdtd.7.3**. The deferral is not "haven't got around to
it"; it is "energy balance has to be re-derived" (Chevalier 1997's
closing remarks acknowledge the analogous incompatibility for the
contemporary Berenger split-field PML).

## 9. TF/SF interaction

Yee already ships a region-coupling mechanism with similar interface
mathematics: the **total-field / scattered-field** decomposition in
Phase 2.fdtd.5 (see [`fdtd-details.md`](./fdtd-details.md) §4).
TF/SF separates a total-field interior from a scattered-field
exterior across a closed rectangular box; the interface is six
axis-aligned faces; an analytic incident-field correction is added
to the Yee curl term on each face. The analogy is direct:

|                       | TF/SF                           | Subgridding                  |
|-----------------------|---------------------------------|------------------------------|
| Two regions           | total / scattered field         | coarse parent / fine nest    |
| Interface             | rectangular box, 6 faces        | rectangular box, 6 faces     |
| Coupling field        | analytic incident E, H          | numerical neighbouring grid  |
| Coupling type         | additive correction per face    | Dirichlet-style overwrite    |
| Reciprocity guarantor | analytic incident exact         | discrete transpose `T = R^T` |

`SubgriddedSolver::step`'s seven-stage sequence is the direct
counterpart of the TF/SF correction sequence in
`crates/yee-fdtd/src/tfsf.rs`. The implementation plan's API sketch
mirrors the TF/SF crate: configuration struct (`SubgridRegion` ↔
`TfsfBox`), solver wrapper (`SubgriddedSolver` ↔ `TfsfDriver`),
per-stage update inside the parent solver's step body. TF/SF and
subgridding cannot share a face in v0 (same auxiliary-state argument
as §8). Phase 2.fdtd.7.0 forbids overlap; Phase 2.fdtd.7.3 takes it
on.

## 10. Validation: the fdtd-007 dielectric-loaded thin-slot case

The walking-skeleton production gate is **fdtd-007**: a thin slot
antenna in an infinite PEC ground plane backed by a dielectric slab.
Geometry (Maloney & Smith 1993): slot `w = 0.5 mm` × `L = 30 mm`;
dielectric `ε_r = 2.2`, thickness `h = 1.524 mm`; delta-gap voltage
drive at the slot centre. Coarse `dx_coarse = 1 mm`; fine nest
`dx_fine = 0.5 mm` over a `(40 × 6 × 4) mm` box centred on slot and
substrate.

**Why this case validates subgridding specifically.** The slot width
`w = 0.5 mm` is exactly **one fine cell** — half a coarse cell on a
uniform `dx = 1 mm` grid, which staircases to zero or one cell and
corrupts the resonant behaviour. A globally uniform `dx = 0.5 mm`
grid resolves the slot correctly but pays `2³ · 2 = 16×` the cost
relative to the coarse-plus-nest split. The fine region covers slot
and immediate substrate; surrounding air, carrying radiating
wavefronts but no sub-wavelength features, stays coarse — the
canonical subgridding sweet spot.

**Reference.** Maloney & Smith, IEEE T-AP 41.5 (1993), pp. 668–676,
Fig. 9. Tolerance: resonant frequency within ±2%, `|S_11|` at
resonance within ±1 dB. A second sanity check runs the same problem
on a **globally uniform `dx = 0.5 mm` grid** and requires agreement
within 0.3% / 0.3 dB. Run-time budget: `< 30 min` `--release` on the
CI Linux runner; if overrun, hardware-gated behind `#[ignore]` per
the mom-001 / Phase 1.5 precedent.

## 11. Limitations and roadmap

The walking-skeleton scope is intentionally narrow; each restriction
maps to a deferred sub-phase.

- **Refinement ratios `r ≥ 3`** (Phase 2.fdtd.7.4). The
  energy-balance transpose for `r = 3, 5, 7` requires a more careful
  weight derivation; weight tables grow with `r²` per face. Higher
  ratios also expose second-order Yee dispersion at the interface.
- **Multi-region nesting and nest-inside-nest** (Phase 2.fdtd.7.1).
  Disjoint nests are independent bookkeeping; nest-inside-nest is
  more delicate (the inner nest sees its outer parent's interface
  as an interpolated source).
- **Dispersive ADE inside the fine region** (Phase 2.fdtd.7.2).
  Auxiliary polarisation variables
  ([`fdtd-details.md`](./fdtd-details.md) §5) have their own
  time-step and re-staggering across the interface, analogous to
  the CPML problem in §8.
- **CPML and TF/SF co-location** (Phase 2.fdtd.7.3). Forbidden in
  v0; the next sub-phase designs the auxiliary-state transpose
  operators.
- **GPU parallelism across the coarse / fine boundary.** Single-GPU
  is straightforward. Multi-GPU (Phase 4) where the fine region
  straddles a slab-decomposition boundary is harder — the
  seven-stage interleaving must align with the halo schedule.

## 12. References

- Okoniewski, M., Okoniewska, E., and Stuchly, M. A. "Three-
  dimensional subgridding algorithm for FDTD." *IEEE Trans. Antennas
  Propag.* 45.3 (1997), pp. 422–429. — Foundational temporal
  sub-cycling pattern; §6 follows it directly.
- Chevalier, M. W., Luebbers, R. J., and Cable, V. P. "FDTD local
  grid with material traverse." *IEEE Trans. Antennas Propag.* 45.3
  (1997), pp. 411–421. — Spatial-interpolation prescription and
  energy-balance closure (the `T = R^T` argument in §5).
- Berenger, J.-P. "Reflection from Courant–Friedrichs–Lewy condition
  violation when a wave penetrates the interface between two FDTD
  grids." *IEEE Trans. Antennas Propag.* 51.10 (2003),
  pp. 2884–2893. — Canonical analysis of the late-time
  asymmetric-coupling instability §5 cancels.
- Berenger, J.-P. "A Huygens subgridding for the FDTD method."
  *IEEE Trans. Antennas Propag.* 54.12 (2006), pp. 3797–3804. —
  Stronger-stability fallback; v0 does not ship it.
- Maloney, J. G., and Smith, G. S. "A study of transient radiation
  from the Wu-King resistive monopole — FDTD analysis and experimental
  measurements." *IEEE Trans. Antennas Propag.* 41.5 (1993),
  pp. 668–676. — fdtd-007 gate; Fig. 9 supplies the reference curve.
- Taflove, A., and Hagness, S. C. *Computational Electrodynamics.*
  3rd ed., Artech House, 2005. §11 — subgridding survey.
