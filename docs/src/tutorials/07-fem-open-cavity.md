# Open-boundary FEM driven sweep from Python

This tutorial walks through Yee's **Phase 4.fem.eig.2** open-boundary
FEM driver from Python. It is the open-region follow-on to
[Tutorial 6 — Dispersive FEM cavity eigenmode](06-fem-cavity-dispersive.md):
where Phases 4.fem.eig.0/1 returned cavity eigenpairs on a fully
PEC-bounded box, the v2 solver consumes a mesh whose exterior faces
are partitioned into PEC, **1st-order Engquist–Majda absorbing
boundary (ABC)**, and **modal wave-port** classes, assembles a
complex-symmetric driven system per swept frequency, and returns the
frequency-swept **S-parameter matrix** via a single complex sparse LU
back-substitution per frequency point.

## Goal

Solve the driven open-region Helmholtz problem

```text
∇ × (1/μ_r) ∇ × E − k₀² ε_r E = 0    in Ω
n̂ × E = 0                            on Γ_PEC
n̂ × ∇×E = −j k₀ n̂ × (n̂ × E)         on Γ_ABC
E_t = (a_inc + b) · e_mode(x,y)      on Γ_port
```

on a small WR-90 stub, recover `|S_{11}(f)|` across an 8–12 GHz sweep,
and confirm the passivity + monotonicity sanity gates baked into the
`fem-eig-003` validation driver. The Phase 4.fem.eig.2 spec
([`docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`](../../superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md))
scopes the v2 deliverable to:

- **1st-order Engquist–Majda ABC** on tagged exterior faces — adds
  `+ j k₀ · area · (n̂ × N_i) · (n̂ × N_j)` per face into the global
  complex stiffness matrix `K(ω)`, promoting the system to
  complex-symmetric even with real `ε_r`.
- **Single dominant-mode wave-port faces** sourced from
  `NumericalCrossSection::e_tangential_at` (Phase 1.3.1.1) with an
  incident modal amplitude on the right-hand side and `S_{11}(f)`
  extracted via modal projection of the FEM solution.
- **Frequency sweep** via per-frequency complex sparse LU
  back-substitution on the Phase 4.fem.eig.1 `faer` surface.
- **Validation gate `fem-eig-003` (WR-90 stub + ABC)**: passivity
  `|S_{11}| ≤ 1 + ε_num` and adjacent-bin smoothness, both enforced
  in default CI; the strict `[−45, −35] dB` ABC absorption window is
  `#[ignore]`'d at the v0 walking-skeleton mesh resolution pending a
  modal-RHS scaling fix (queued as Phase 4.fem.eig.2.0.1 /
  4.fem.eig.2.5).

The v2 scope, deferrals, and the absorption-floor finding are
recorded in
[ADR-0040](../decisions/0040-phase-4-fem-eig-2-open-boundary-scope.md);
the parallel pattern of open-question surfacing (rather than ad-hoc
gate weakening) is also tracked in
[ADR-0041](../decisions/0041-fdtd-007-reference-correction.md).

## Prerequisites

- Rust 1.92+ (`rust-toolchain.toml` pins the toolchain; the Python
  wheel is built from source).
- Python 3.10 through 3.14 — the wheel is `abi3-py310`, so any
  interpreter in that range works without rebuilding.
- `uv pip install maturin numpy pytest` — or the equivalent
  `pip install` invocation.

Install matches [Tutorial 6](06-fem-cavity-dispersive.md) exactly:

```bash
uv venv .venv
source .venv/bin/activate
uv pip install maturin numpy pytest
cargo build --release -p yee-fem -p yee-py
cd crates/yee-py
maturin develop --release
python -c "import yee.fem; help(yee.fem.solve_open_cavity)"
```

`maturin develop --release` is non-negotiable on the open-boundary
path: each swept frequency point re-assembles a complex sparse
pencil and factors a `Complex64` LU. Debug builds are an order of
magnitude slower and unsuitable for sweep work.

## WR-90 stub — first call

The minimal end-to-end snippet — air-filled WR-90 stub, ABC at
`z = 0`, TE_{10} wave-port at `z = d`, sweep across 8–12 GHz:

```python
import math
import yee.fem

a = 0.02286    # WR-90 broad wall (m)
b = 0.01016    # WR-90 narrow wall (m)
d = 0.030      # axial stub length (m)

# Constant tangential TE_{10} amplitude at the broad-wall midpoint.
# e_mode(x,y) = ŷ · sqrt(2 / (a·b)) · sin(π x / a); sampled at x = a/2
# the sin term is 1, so modal_e_t reduces to (0, sqrt(2/(a·b)), 0).
te10_norm = math.sqrt(2.0 / (a * b))

# 50-point uniform sweep, 8-12 GHz.
omegas = [(8.0 + 4.0 * i / 49.0) * 1e9 for i in range(50)]

s = yee.fem.solve_open_cavity(
    a, b, d, 8, 4, 16,
    materials=[{"tag": 0, "eps_inf": 1.0, "mu_r": 1.0, "poles": []}],
    port_faces=[
        {"axis": "z", "side": "high", "port_id": 0,
         "modal_e_t": (0.0, te10_norm, 0.0)},
    ],
    abc_faces=[{"axis": "z", "side": "low"}],
    omegas_hz=omegas,
)

for i, om in enumerate(omegas):
    s11_db = 20 * math.log10(abs(s[i, 0, 0]))
    print(f"f = {om/1e9:.2f} GHz   |S11| = {s11_db:+.3f} dB")
```

`s` is a complex `numpy.ndarray` of shape `(n_omegas, n_ports,
n_ports)` — here `(50, 1, 1)` — indexed as `s[i, k, j] = S_{kj}(ω_i)`.

Expected output: `|S_{11}(f)| ≈ 0 dB` (numerical-1.0) across the
entire band. That is **not** the published 1st-order Engquist–Majda
floor (`~ −40 dB`) — it is the v0 walking-skeleton modal-RHS scaling
saturation surfaced by Track BBBBBBBBB during the `fem-eig-003`
landing (see "Known limitations" below). The shape of the sweep is
flat-at-unity rather than noisy, which is itself a useful diagnostic:
the structure is being treated as effectively closed, not as a leaky
absorber.

## Face classification

The Phase 4.fem.eig.2 binding partitions every exterior face into one
of three tags. Unannotated exterior faces default to **PEC** — every
Phase 4.fem.eig.0/1 caller therefore round-trips unchanged.

- **`port_faces`** — one dict per physical wave-port face. Required
  keys:
  - `axis`: `"x" | "y" | "z"` — the axis-aligned face normal direction.
  - `side`: `"low" | "high"` — which of the two opposing faces on
    that axis.
  - `port_id`: `int` — the column index in the returned S-parameter
    matrix.
  - `modal_e_t`: `(ex, ey, ez)` — constant tangential modal `E`-field
    sampled at the face centroid. For TE_{10} on a rectangular
    waveguide the sampled value at `x = a/2` is the orthonormalised
    amplitude `sqrt(2 / (a·b))` along the broad-wall-perpendicular
    axis.
- **`abc_faces`** — one dict per Engquist–Majda absorbing face. Keys
  `axis` and `side` only; the ABC face block is parameter-free at v0
  (the `+ j k₀` weight is implied).
- **Default-PEC** — every exterior face not listed under `port_faces`
  or `abc_faces` is treated as a tangential-`E`-zero Dirichlet wall.
  This is the same boundary condition the Phase 4.fem.eig.0
  closed-cavity solver enforces on every exterior face.

> **Multi-port S-matrix extraction is deferred.** v0 supports any
> number of port faces in the `port_faces` list but the returned
> S-parameter columns beyond `port_id = 0` are populated from a
> *single-incident* solve per frequency, not the per-port excitation
> matrix the multi-port S-formalism requires. Full multi-port driving
> lands in **Phase 4.fem.eig.2.0.2**; until then, treat `n_ports > 1`
> sweeps as `S_{*0}(f)` columns only.

## Validation

The corresponding validation gate, **`fem-eig-003`**, lives at
`crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` and drives
the public `yee_validation::run_fem_eig_003_wr90_stub_abc` helper.
The driver assembles the spec-scale mesh `(nx, ny, nz) = (16, 8, 24)`
(18 432 tets) with the analytic TE_{10} modal source on the `z = d`
face and the 1st-order ABC on `z = 0`, then sweeps 50 uniform points
across 8–12 GHz.

The default-CI gates are:

1. **Passivity** — `|S_{11}(f)| ≤ 1 + ε_num` with `ε_num = 0.05` at
   every swept point (the strict `< 1` continuum bound is `#[ignore]`'d
   under the plan E5 escape hatch — see below).
2. **Smoothness** — adjacent-bin `|Δ(20·log₁₀|S_{11}|)| ≤ 10 dB`. No
   spurious near-cavity resonance appears in the propagating band.
3. **Finiteness** — every entry is finite (no NaN / Inf).

The strict ABC-absorption gate `|S_{11}(f)|_dB ∈ [−45, −35] dB` is
`#[ignore]`'d pending the Track CCCCCCCCC modal-RHS scaling fix.
Run it explicitly via `cargo test -- --ignored` once that lands.

The row in `crates/yee-fem/validation/README.md` carries the
spec-scale mesh, the per-gate tolerance breakdown, and the
finding-history pointers.

## CFS-PML mode (Phase 4.fem.eig.3.5)

Phase 4.fem.eig.3.5 ships **CFS-PML** (Complex Frequency Shifted
Perfectly Matched Layer, Roden-Gedney 2000) as an alternative
truncation kernel to the Engquist-Majda surface integral. CFS-PML
is a thin (default 6-cell) volumetric buffer of additional tetrahedra
outside the original cavity, in which the constitutive tensor becomes
the stretched-coordinate form `ε(ω) · Λ(ω)`. Unlike the local
surface-integral ABC, the volumetric PML absorbs **off-normal** and
**evanescent** modal content uniformly — exactly the regime where
the 2nd-order Engquist-Majda operator's intrinsic floor exposes
itself on the WR-90 stub (`|S_{11}|` saturates near `-2e-2 dB` at the
`(24, 12, 36)` mesh tier; see the fem-eig-003 strict-gate test
docstring for the ADR-0042 §risks deferral path).

To engage CFS-PML from the Python binding, pass a `pml_config` dict
to `yee.fem.solve_open_cavity`:

```python
import yee.fem

a, b, d = 0.02286, 0.01016, 0.030  # WR-90 stub
nx, ny, nz = 24, 12, 36

s = yee.fem.solve_open_cavity(
    a, b, d, nx, ny, nz,
    materials=[],
    port_faces=[{
        "axis": "z",
        "side": "high",
        "port_id": 0,
        "modal_e_t": (0.0, 1.0, 0.0),
    }],
    abc_faces=[{"axis": "z", "side": "low"}],
    omegas_hz=[10.0e9],
    coupled_whitney=True,
    pml_config={
        "thickness_cells": 6,   # default; 6 to 10 per Roden-Gedney 2000 §III
        "kappa_max": 5.0,       # coordinate stretching at outer surface
        "m": 3,                 # polynomial grading order
        # sigma_max and alpha_max default to 0.0 → auto-resolved
    },
)
```

The `pml_config` keys mirror `yee_fem::PmlConfig`. When supplied,
`pml_config` supersedes `abc_order` — the volumetric PML absorbs in
the bulk and the surface-integral kernel is suppressed on every
`abc_faces` entry. Cartesian-aligned PML only in v3.5 (ADR-0043 §4);
multi-axis edge / corner wedges land in Phase 4.fem.eig.3.5.1.

### Grading parameter defaults (Phase 4.fem.eig.3.5.2)

The v3.5.2 H4 ablation sweep (33 configurations: H1 baseline + H2
κ_max sweep + H3 (m, thickness) sweep + H4 α-grading sweep) retired
the fem-eig-003 strict absorption-floor gate. Winning defaults:

| Field | Value | Why |
|-------|-------|-----|
| `thickness_cells` | **16** | extended H3 ablation (v3.5.1 capped at 10) |
| `m` | **3** | polynomial-grading sweet spot for 16-cell shell |
| `kappa_max` | **2** | smaller than Roden-Gedney 2000 Table I default 5 (per Berenger 2002 §V) |
| `alpha_grading_order` | **1** | linear `α_α(d) = α_max·(1-d/D)` per Berenger 2002 §VI |
| `sigma_max`, `alpha_max` | `0.0` (auto-resolved) | per-axis from `PmlMeshMeta` (v3.5.1 R1) |

Measurement on fem-eig-003 WR-90 stub: `|S_{11}(f)|` band `[-71.53,
-55.58] dB` — ~50 dB total improvement over v3.5 baseline `[-11.0,
-7.48] dB`. Both `fem_eig_003_strict_absorption_floor_gate` and
`fem_eig_003_strict_passive_bound_continuum_limit` un-ignore and
pass in default CI.

#### Knob → effect table (full H4 sweep, fem-eig-003 WR-90 stub)

| Hypothesis | `κ_max` | `m` | `t` | `α_order` | `s11_db_min` | `s11_db_max` |
|------------|---------|-----|-----|-----------|--------------|--------------|
| H1 (per-axis baseline)      | 5 | 3 | 6  | 0 | -31.20 | -21.74 |
| H2 κ_max sweep              | 2 | 3 | 6  | 0 | -30.90 | -22.55 |
| H3 most-aggressive probe    | 2 | 4 | 10 | 0 | -58.13 | -35.45 |
| **H4 winner (default)**     | **2** | **3** | **16** | **1** | **-71.53** | **-55.58** |
| H4 best per m=4             | 2 | 4 | 16 | 0 | -82.73 | -54.23 |

H4 winner is `(κ=2, m=3, t=16, α_order=1)` — the smallest
`m × t × (α_order + 1)` product retiring fem-eig-003 with comfortable
margin past the -40 dB upper bound.

#### fem-eig-006 status

The 100:10:1 high-aspect cavity at 30 GHz remained at `|S_{11}| =
0.926` across all 18 H4 rows. α grading is orthogonal to that
fixture — dominant modal content is not normal-incidence at the +x
face. `fem_eig_006_magnitude_bounded` stays `#[ignore]`'d; queued
for Phase 4.fem.eig.3.5.3 / 4.fem.eig.4 (rotated PML / wave-port
termination / multi-face wedges).

### Phase 4.fem.eig.3.5.3 — fem-eig-006 wave-port retirement attempt (W1)

ADR-0046 (`docs/src/decisions/0046-phase-4-fem-eig-3-5-3-fem-eig-006-retire.md`)
chose **W1**: replace the +x CFS-PML shell with a TE_{10}
`FaceKind::WavePort(1)` termination per Jin §10.6. Berenger 1996
§IV-A predicted Cartesian-aligned CFS-PML cannot absorb guide-modes
at the +x face — the 18-row H4 sweep (frozen at `|S_{11}| = 0.926`)
confirmed it.

The W1 swap is a ~10-line driver edit
(`run_fem_eig_006_high_aspect_pml` in `crates/yee-validation/src/lib.rs`)
that drops `extend_mesh_with_pml` and the CFS-PML rim altogether,
running the native `(16, 3, 2)` Kuhn-6 cavity (576 tets vs the
v3.5.2 ~580 extended tets).

**Measurement on W1**:

| Configuration                | `|S_{11}|(30 GHz)` | dB     |
|------------------------------|--------------------|--------|
| v3.5.2 CFS-PML (best H4 row) | 0.926              | -0.67  |
| v3.5.3 W1 TE_{10} wave-port  | 0.925644           | -0.67  |

The W1 termination matches the PML measurement within numerical
noise. Per ADR-0046 §Decision (5) and spec §7 (a): the broad-wall
`b = 10 mm` puts the **TE_{20} cutoff at exactly 30 GHz**, so a
TE_{10}-only port cannot project out the dominant higher-order
modal content driving reflection.

`fem_eig_006_magnitude_bounded` therefore stays `#[ignore]`'d;
tolerance `< 0.1` is **not** weakened. Multi-mode wave-port
extension (add `TE_{20}` and `TE_{01}` mode shapes to the +x
`PortDefinition`) is queued for **Phase 4.fem.eig.3.5.4**.
fem-eig-003 (CFS-PML absorption floor on the cubic cavity) remains
retired at the v3.5.2 strict band `[-71.53, -55.58] dB` — that
driver is untouched.

### Phase 4.fem.eig.3.5.4 — multi-mode wave-port attempt + cutoff-degeneracy finding

ADR-0047 (`docs/src/decisions/0047-phase-4-fem-eig-3-5-4-multi-mode-wave-port.md`)
extended `PortDefinition` from a single
`(beta_mode, modal_e_t)` closure pair to a modal basis
`{ modes: Vec<PortMode> }`, with the `single_mode` constructor
preserving v3.5.3 call shape. The v3.5.4 driver populates the +x
face with three modes (`TE_{10}` driving, `TE_{20}` and `TE_{01}`
outgoing-only):

| Configuration                                  | `|S_{11}|`  | dB     |
|------------------------------------------------|-------------|--------|
| v3.5.2 CFS-PML (best H4 row)                   | 0.926       | -0.67  |
| v3.5.3 W1 single-mode TE_{10}                  | 0.925644    | -0.67  |
| v3.5.4 multi-mode (TE_{10}, TE_{20}, TE_{01})  | 0.925637    | -0.67  |

The multi-mode basis collapses to single-mode at 30 GHz exactly
because **`TE_{20}` cutoff sits at the test frequency**:

* Port-face broad wall is `B = 10 mm` (the cavity's 10 : 1 : 1
  narrow direction at the port face, not the propagation length
  `A = 100 mm` — the v3.5.4 design spec §2.2 mis-derived this).
* `TE_{20}` cutoff `f_c = c / B = c / 0.010 m = 30.0 GHz` exactly.
  At cutoff `β = 0`, the per-mode stiffness block contribution
  vanishes identically.
* `TE_{01}` cutoff `f_c = c / (2 D) = c / 0.002 m = 150.0 GHz` —
  evanescent at 30 GHz; carries no propagating modal content.

The multi-mode infrastructure (`PortMode`, `Vec<PortMode>`,
multi-mode summation in `scatter_port_face`, driving-mode
selection in S-parameter extraction) lands cleanly and will earn
its keep at any non-degenerate test frequency. The
`fem_eig_006_magnitude_bounded` gate stays `#[ignore]`'d and
tolerance `< 0.1` is **not** weakened. ADR-0048 (Phase
4.fem.eig.3.5.5) queues two candidate dispositions:

(a) **Retune `FEM_EIG_006_F_HZ`** to a frequency well above the
    `TE_{20}` cutoff (e.g. 40 GHz) so the multi-mode basis carries
    real propagating modal content and the projection step has
    something non-trivial to do.

(b) **Land an absorbing-mode wave-port** (Lee-Mittra 1997) that
    handles arbitrary modal content via evanescent absorption
    rather than modal projection.

### Phase 4.fem.eig.3.5.5 — frequency retune to 40 GHz (ADR-0048 Option (a))

Phase 4.fem.eig.3.5.5 executed ADR-0048 Option (a): `FEM_EIG_006_F_HZ`
was retuned from 30 GHz to **40 GHz**, where `TE_{20}` propagates with
`β ≈ 554 rad/m` (33% above its `f_c = c / B = 30 GHz` cutoff) instead of
sitting exactly on it. The cavity, mesh, and 3-mode driver are otherwise
unchanged, so the multi-mode basis now carries real propagating content
for the projection step to terminate:

| Configuration                                  | `|S_{11}|`  | dB     |
|------------------------------------------------|-------------|--------|
| v3.5.2 CFS-PML (best H4 row)                   | 0.926       | -0.67  |
| v3.5.3 W1 single-mode TE_{10}                  | 0.925644    | -0.67  |
| v3.5.4 multi-mode (TE_{10}, TE_{20}, TE_{01}) @ 30 GHz | 0.925637 | -0.67  |
| **v3.5.5 multi-mode @ 40 GHz (TE_{20} propagating)** | **0.955397** | **-0.40** |

The retune **did not retire** the gate — `|S_{11}|(40 GHz) = 0.955397`
is marginally *above* the cutoff-degenerate 30 GHz value, so the v3.5.4
modal degeneracy was not the binding constraint. A one-shot refinement
probe (transverse mesh `NY 3→9`, `NZ 2→6`, 5184 tets) gave
`|S_{11}| = 0.913956` (-0.78 dB): a 9× transverse element-count increase
moved the residual only ~0.04, so the reflection is **not**
discretisation-limited either (probe reverted, native (16, 3, 2) mesh
stands).

**Disposition: escape-hatch.** With both modal degeneracy (v3.5.4) and
discretisation (this probe) excluded, the residual ~0.95 is a genuine
limitation of the modal-projection wave-port — projecting onto a finite
`TE_{mn}` basis cannot fully match the field at the truncation face of a
strongly off-square (100 : 10 : 1) cavity. The `fem_eig_006_magnitude_bounded`
gate stays `#[ignore]`'d (tolerance `< 0.1` **not** weakened), and ADR-0049
(`docs/src/decisions/0049-phase-4-fem-eig-3-5-6-absorbing-mode-wave-port.md`)
queues **Phase 4.fem.eig.3.5.6** to land the Lee-Mittra 1997 absorbing-mode
wave-port (ADR-0048 Option (b)). `FEM_EIG_006_F_HZ = 40 GHz` is retained as
the operating point for that work — reverting to 30 GHz would reintroduce
the cutoff degeneracy.

#### Per-fixture override pattern

Users who want to revert to the v3.5.0 OOOOOOOOO baseline (thinner
PML, no α grading) can override the default:

```python
import yee.fem

s = yee.fem.solve_open_cavity(
    a, b, d, nx, ny, nz,
    materials=[],
    port_faces=[{"axis": "x", "side": "high", "port_id": 0,
                 "modal_e_t": (0.0, 1.0, 0.0)}],
    abc_faces=[{"axis": "x", "side": "low"}],
    omegas_hz=[30.0e9],
    coupled_whitney=True,
    pml_config={
        "thickness_cells": 6,         # back to OOOOOOOOO default
        "kappa_max": 5.0,
        "m": 3,
        "alpha_grading_order": 0,     # constant α_max (v3.5.1)
    },
)
```

## Known limitations

- **`|S_{11}(f)| = 1.0` saturation at v0 mesh resolution.** Track
  BBBBBBBBB's E5 landing measured `|S_{11}(f)| ≈ 1.000_000_000`
  numerically across the entire 8–12 GHz sweep on the spec-scale
  `(16, 8, 24)` mesh — the published 1st-order Engquist–Majda
  reflection floor `~ −40 dB` does **not** resolve at the
  walking-skeleton mesh + face-centroid modal-RHS combination shipped
  in v0. The ABC face block measurably differs from PEC (`Im(S_{11})`
  differs at `~1e-10` vs `~1e-8` on the coarse mesh), but the
  modal-source RHS is too weak to discriminate from a fully-PEC
  structure on a real WR-90 sweep. The strict absorption-floor gate
  is `#[ignore]`'d under the plan E5 escape hatch ("if walking-skeleton
  physics doesn't resolve `-40 dB` at 25 k tets, document and
  continue"). Full discussion in `crates/yee-fem/validation/README.md`
  E5 findings.
- **Track CCCCCCCCC modal-RHS scaling fix is in flight.** The agreed
  follow-up is to re-derive the per-Gauss-point modal projection
  (cubic interpolation per ADR-0040 §C-3 amendment) and re-scale the
  RHS to match the spec §4.3 closed form. Once it lands, the
  `#[ignore]`'d strict gate can be lifted with a single attribute
  removal — no API surface change is required, and every snippet in
  this tutorial keeps producing the same `s` array shape.
- **Single dominant mode per port at v0.** Higher-order modes are
  *captured* in the modal reflection spectrum (they show up as
  band-edge artefacts past the TE_{10} cutoff at ~6.56 GHz on WR-90)
  but are not *driven*. Multi-mode incident excitation lands in
  4.fem.eig.2.0.2.
- **Combining the v2 open-boundary path with v1's dispersive Newton
  tracker** is a Phase 4.fem.eig.2.1 superposition exercise. The
  driven sweep at v0 supports lossless or constant-real-loss interior
  media only; do not mix `poles=[…]` materials with `port_faces` in
  the same call (the call accepts the shape, but the dispersive
  ε(ω) is not iterated and the result is not physically meaningful).

## What's next

Phase 4.fem.eig.2 is a strict extension of the closed-cavity walking
skeleton — every Phase 4.fem.eig.0 / 4.fem.eig.1 caller round-trips
unchanged, and the `OpenBoundarySolver` is rejected at construction
with an empty `abc_faces` and empty `ports` (no excitation → no
well-posed driven problem). The roadmap lays out the immediate
follow-ons:

- **Phase 4.fem.eig.2.0.1** — cubic / per-Gauss-point modal-profile
  interpolation on the FEM port face, retiring the v0 face-centroid
  sample and lifting the strict `fem-eig-003` absorption gate.
- **Phase 4.fem.eig.2.0.2** — multi-port incident-excitation matrix
  formalism; full `n_ports × n_ports` S-matrix per swept frequency.
- **Phase 4.fem.eig.2.5** — 2nd-order Engquist–Majda / Higdon /
  CFS-PML termination if the 1st-order floor cannot meet a future
  published-benchmark tolerance.
- **Phase 4.fem.eig.3** — driven-Maxwell solver: dielectric-resonator
  antennas, coax-fed dipoles in ABC-terminated FEM boxes
  (`fem-eig-004` cross-checks against the `mom-001` NEC-4 reference),
  iris-coupled bandpass filters.

If you want the full design rationale — surface-term derivation,
modal projection convention, port-vs-PEC edge precedence, ABC face
orientation — read the spec at
[`docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`](../../superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md)
and the plan at
[`docs/superpowers/plans/2026-05-19-phase-4-fem-eig-2-open-boundary.md`](../../superpowers/plans/2026-05-19-phase-4-fem-eig-2-open-boundary.md).

## References

- **Engquist, B. and Majda, A.**, "Absorbing boundary conditions for
  the numerical simulation of waves", *Math. Comp.* 31 (1977),
  pp. 629–651 — the canonical 1st-order ABC derivation that v0 ships.
- **Jin, J.-M.**, *The Finite Element Method in Electromagnetics*,
  3rd ed., Wiley 2014 — Ch. 10 (driven FEM analysis), §10.4 (ABC
  face contributions), §10.5 (wave-port modal decomposition), §10.7
  (S-parameter extraction).
- **Pozar, D. M.**, *Microwave Engineering*, 4th ed., Wiley 2012 —
  §3.3 (waveguide TE/TM modes, propagation constants, closed-form
  modal characterisation of a uniformly-terminated rectangular
  waveguide section in the dominant-mode band).
- **ADR-0040** —
  [`docs/src/decisions/0040-phase-4-fem-eig-2-open-boundary-scope.md`](../decisions/0040-phase-4-fem-eig-2-open-boundary-scope.md):
  v2 scope (1st-order Engquist–Majda + single-mode modal port), PML /
  2nd-order ABC deferral, `fem-eig-003` `[−45, −35] dB` window and
  the absorption-floor finding.
- **ADR-0041** —
  [`docs/src/decisions/0041-fdtd-007-reference-correction.md`](../decisions/0041-fdtd-007-reference-correction.md):
  the parallel "surface the open question as an ADR rather than
  weaken the gate" pattern this tutorial's `#[ignore]`'d strict gate
  also follows.
- **Yee project spec** —
  [`docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`](../../superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md).
- **Yee project plan** —
  [`docs/superpowers/plans/2026-05-19-phase-4-fem-eig-2-open-boundary.md`](../../superpowers/plans/2026-05-19-phase-4-fem-eig-2-open-boundary.md);
  E1 (ABC face block), E2 (modal RHS), E3 (driven assembly), E4
  (sweep + S-parameter extraction), E5 (`fem-eig-003` gate), E6
  (`yee.fem.solve_open_cavity` Python binding).
- **Phase 1.3.1.1** —
  `crates/yee-mom/src/eigensolver/`: the 2-D Nedelec cross-section
  eigensolver whose `NumericalCrossSection::e_tangential_at` accessor
  is the FEM port's modal-profile source.
- **Phase 4.fem.eig.0** —
  [Tutorial 4 — FEM cavity eigenmode from Python](04-fem-cavity-eigenmode.md);
  the lossless closed-cavity walking skeleton this tutorial extends.
- **Phase 4.fem.eig.1** —
  [Tutorial 6 — Dispersive FEM cavity eigenmode from Python](06-fem-cavity-dispersive.md);
  the lossy-dispersive sibling whose complex sparse LU surface v2 reuses.
- **Validation driver** —
  `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`
  (Track BBBBBBBBB landing, E5).
- **Python test fixture** —
  `crates/yee-py/tests/test_fem_open_boundary.py` — the pytest
  re-running the smoke + passivity + ABC-monotonicity gates from
  Python.
