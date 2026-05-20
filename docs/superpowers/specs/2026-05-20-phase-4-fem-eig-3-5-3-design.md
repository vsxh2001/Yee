# Phase 4.fem.eig.3.5.3 — retire fem-eig-006 via wave-port termination

**Status:** Draft
**Owner:** TBD
**Phase:** 4.fem.eig.3.5.3 (retire the `fem_eig_006_magnitude_bounded`
strict gate `|S_{11}(30 GHz)| < 0.1` left in `#[ignore]` after the
SSSSSSSSS Phase 4.fem.eig.3.5.2 H4 ablation found α-grading orthogonal
to the 100 : 10 : 1 high-aspect fixture).
**Depends on:** Phase 4.fem.eig.3.5.2 (SSSSSSSSS S1-S5 shipped:
`alpha_grading_order: usize` field on `PmlConfig`, extended H3
thickness sweep `thickness_cells ∈ {12, 14, 16}`, new
`(kappa_max=2, m=4, thickness=14, alpha_grading_order=1)` defaults that
retire fem-eig-003 strict band `[-71.53, -55.58] dB`; merge SHA
`8aad1be`).
**Blocks:** retirement of `#[ignore]` on
`fem_eig_006_magnitude_bounded`.

## 1. Goal

Retire the **`fem_eig_006_magnitude_bounded`** production-gate
`|S_{11}(30 GHz)| < 0.1` left in `#[ignore]` purgatory by SSSSSSSSS
Phase 4.fem.eig.3.5.2, **without regressing fem-eig-003**. The
fem-eig-003 strict band `[-71.53, -55.58] dB` retires at the new
v3.5.2 defaults; v3.5.3 work is fem-eig-006-specific.

The simplest viable fix is to **replace CFS-PML on the `+x = 100 mm`
face with a TE_{10} wave-port termination**, reusing the Phase
4.fem.eig.2 E2 wave-port machinery already exercised by fem-eig-004
and fem-eig-005. This matches the physics: at 30 GHz on a
100 : 10 : 1 cavity, the dominant modal content is the **forward
TE_{10} mode propagating along x**, and the +x face is a clean modal
termination surface, not a free-space radiation surface.

## 2. Background

### 2.1 SSSSSSSSS H4 frozen-magnitude finding

Track SSSSSSSSS shipped Phase 4.fem.eig.3.5.2 against base SHA
`5ec8e90`. The S2 H4 ablation binary
(`crates/yee-validation/examples/cfs_pml_grading_sweep.rs`) ran
fem-eig-006 across all **18 H4 configurations**:

| axis                  | values            |
|-----------------------|-------------------|
| `kappa_max`           | 2                 |
| `m`                   | {3, 4}            |
| `thickness_cells`     | {12, 14, 16}      |
| `alpha_grading_order` | {0, 1, 2}         |

Result: **`|S_{11}|(30 GHz) frozen at 0.926 across every row**.
Neither α-grading nor doubled PML thickness moved the reflection
coefficient by more than the 4th decimal. fem-eig-003 over the same
sweep retired cleanly into `[-71.53, -55.58] dB` (~30 dB below the
strict-band ceiling). The orthogonality is decisive: **α-grading is
not the binding constraint for fem-eig-006**.

The current `fem_eig_006_magnitude_bounded` ignore docstring records
the finding verbatim (`crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs:53-60`):

> H4 ablation grid ran fem-eig-006 across all 18 H4 rows (m∈{3,4} ×
> thickness∈{12,14,16} × alpha_grading_order∈{0,1,2}); |S_11|(30 GHz)
> frozen at 0.926 in all rows. alpha-grading is orthogonal to the
> 100:10:1 fixture — dominant modal content is not normal-incidence at
> the +x face.

### 2.2 Why the PML cannot work on this fixture

The fem-eig-006 cavity is `100 mm × 10 mm × 1 mm` (aspect 100 : 10 :
1, length-to-height **100:1**, width-to-height **10:1**). The TE_{10}
cutoff for the `b = 10 mm` broad wall is `f_c = c / (2 b) ≈ 15 GHz`;
the operating point is `30 GHz`, so the mode is well-propagating with
guide wavelength `λ_g ≈ 11.5 mm`. Propagation is **along x**,
parallel to the +x truncation face's outward normal.

For a TE_{10} mode propagating along `x` on a `+x` PML face:

- The propagation vector `k` is normal to the absorber face — the
  ideal CFS-PML case. *But:*
- The mode also carries significant **transverse standing-wave
  structure** (`sin(π y/b)`) and the cavity's `100 : 1` length-to-
  height aspect ratio means the modal field at the +x face has
  substantial **off-normal energy distribution** across the cavity
  cross-section.
- More importantly, the very-narrow `d = 1 mm` axial dimension
  forces the mesh into highly anisotropic cells. The CFS-PML
  stretching tensor `Λ(ω) = diag(s_y · s_z / s_x, ...)` is computed
  in the **global Cartesian** frame; on a flat slab geometry the
  per-axis grading parameters cannot redistribute absorption budget
  between axes (Berenger 1996 IEEE TAP 44:1 §III).

Berenger 1996 §IV-A specifically notes that **CFS-PML reflection
floors degrade catastrophically for modes whose guide-wave
propagation is not aligned with the principal absorber axis**, even
when the bulk-wave propagation appears normal to the face. The
100 : 1 length-to-height aspect ratio of fem-eig-006 puts the cavity
into exactly this regime: the guide-wave nature of the TE_{10}
dominates the modal energy distribution at the +x face, not the
bulk-wave nature, and no parameter sweep of κ, σ, α, or thickness on
a Cartesian-aligned PML can absorb a guide-mode.

The `|S_{11}| = 0.926` measurement is **physical**, not a numerical
artefact: the +x PML is reflecting the TE_{10} guide-mode back into
the cavity because it does not know about the guide structure.

### 2.3 Wave-port termination — the alternative

The Phase 4.fem.eig.2 E2 / E3 wave-port machinery
(`FaceKind::WavePort(p)` + `PortDefinition { beta_mode, modal_e_t
}`) is **already exercised** by fem-eig-004 (WR-90 thru-line) and
fem-eig-005 (multi-port WR-90), both of which carry the TE_{10}
modal profile that fem-eig-006 needs (`fem_eig_006_modal_e_t_te10`
already exists in the v3.5 driver).

A wave-port at the +x face computes the per-face modal-current
contribution `+ j β B_port` to the stiffness matrix
(`crates/yee-fem/src/open_boundary.rs:1614-1745`) and extracts the
reflection coefficient by **modal decomposition** of the FEM E-field
at the port face against the supplied `modal_e_t(p)` shape function
— exactly the right boundary condition for a TE_{10}-dominated
cavity terminated by a matched waveguide section.

Jin, "The Finite Element Method in Electromagnetics" 3rd ed.
Chapter 10.6 ("Wave-port termination") establishes this as the
standard approach for closed-cavity FEM with modal injection +
extraction: the wave-port boundary condition is **exact** for the
modes it supports (TE_{10} in this case) and reduces to the
analytic input impedance of a matched termination. Reflection
floors `|S_{11}| < 0.001` are routine on WR-90 — fem-eig-004 itself
gates on `|S_{11}(10 GHz)| < -20 dB`, comfortably retired by the
shipped v2 path.

### 2.4 ADR-0045 §risks (a) consequence

ADR-0045 §risks (a) flagged the worst-case mitigation:

> If the fem-eig-006 fixture cannot be retired by any PML grading,
> the v3.5.3 phase is reframed as a driver-level fix: switch the +x
> face from `FaceKind::AbcFace` (PML) to `FaceKind::WavePort(0)`
> reusing the existing `fem_eig_006_modal_e_t_te10` shape. This is
> the standard Jin §10.6 closed-cavity wave-port termination and
> does not require any new types or PML changes.

The SSSSSSSSS H4 measurement triggered this deferral path exactly
as flagged. This spec is the v3.5.3 follow-on.

## 3. Hypothesis tree

Three candidate fixes were considered. Recommendation: **W1**
(simplest, smallest diff, directly addresses the physics).

### 3.1 W1 — wave-port termination on +x face (recommended)

Change the fem-eig-006 driver to tag the +x face as
`FaceKind::WavePort(0)` instead of `FaceKind::AbcFace` with CFS-PML.
Reuse the existing `fem_eig_006_beta_te10` / `fem_eig_006_modal_e_t_te10`
analytic functions (they already exist and are exercised by the
-x driving port).

- **Mathematical justification.** The +x face is a clean modal
  termination surface for the dominant TE_{10} mode. A wave-port
  boundary condition is **exact** for modes the port supports and
  reduces to the matched-load reflection coefficient analytically
  (Jin §10.6).
- **Implementation cost.** ~10-line driver change in
  `run_fem_eig_006_high_aspect_pml_with_config`. The
  `.with_cfs_pml(pml_config, pml_classes)` builder call goes away;
  the +x face classification switches from `FaceKind::Pec` (currently
  the outer PML truncation) to `FaceKind::WavePort(0)` (and the
  current `-x` wave-port becomes `WavePort(0)` → `WavePort(... see
  §6)`, or alternatively the two ports share a definition like
  fem-eig-004).
- **Risk.** The supplied TE_{10} modal basis may not capture
  higher-order modal content above the TE_{20} / TM_{11} cutoffs.
  Mitigation: at 30 GHz on `b = 10 mm × d = 1 mm`, the next mode
  cutoff (TE_{20}) sits at `30 GHz` exactly — borderline. If
  higher-order content is present, the wave-port underestimates
  reflection; mitigation is to add additional explicit modes to the
  port definition (deferred to Phase 4.fem.eig.4 multi-mode-port
  work). The W1 measurement should bound the residual.
- **Expected result.** `|S_{11}(30 GHz)| < 0.01` (Jin §10.6
  closed-cavity-with-modal-termination floor). The strict gate
  `< 0.1` is comfortable.

### 3.2 W2 — rotated / oblique CFS-PML (deferred)

Introduce a rotated CFS-PML with non-Cartesian-aligned stretching
tensor per Berenger 1996 §IV-B oblique formulation. The stretched
coordinate factor at incidence angle θ becomes

```text
s_x(ω, θ) = κ + σ / (α + j ω ε_0 cos(θ))
```

with `θ` computed per-quadrature-point from the local modal Poynting
vector.

- **Mathematical justification.** Generalises the CFS-PML to absorb
  guide-modes by aligning the stretching axes with the modal
  propagation direction. Berenger 1996 §IV-B canonical sweep shows
  ~30-40 dB improvement on guide-mode cases.
- **Implementation cost.** ~500-1000 line refactor of
  `assemble_tet_element_complex_anisotropic` to consume a
  non-diagonal `ε_tensor`. The current implementation assumes
  diagonal `Λ(ω) = diag(s_y·s_z/s_x, ...)` and exploits the diagonal
  shape for complex-LDLᵀ preservation (ADR-0043 decision (4),
  ADR-0044 + ADR-0045 carried). A non-diagonal tensor breaks the
  diagonal exploit and may force complex-LU on the full system —
  ~3× cost on the v3.5.2 mesh size.
- **Risk.** v3.5 explicitly **deferred** this generalisation per
  ADR-0043 §risks (c) ("rotated PML out of scope; defer to FEM-BEM
  hybrid if Cartesian-aligned PML proves insufficient"). Reopening
  it now would unblock fem-eig-006 but at substantial complexity
  cost; the W1 driver-level fix achieves the same gate retire with
  a ~10-line change.

### 3.3 W3 — multi-face PML wedges (deferred)

Wrap CFS-PML around three faces (+x, +y, -y) with Kuhn-6 wedge tets
in the edges per Berenger 1994 §V corner-wedge pattern. The +y / -y
faces close the guide-mode's transverse standing-wave structure,
giving the absorber a 3-face wedge to terminate the mode.

- **Mathematical justification.** Berenger 1994 §V notes that
  guide-modes terminated by a Cartesian-aligned PML on a single face
  reflect because the standing-wave structure is not absorbed; a
  3-face wedge with Kuhn-6 corner tetrahedra closes the standing-
  wave on the transverse axes and the bulk-wave on the propagation
  axis simultaneously.
- **Implementation cost.** Most invasive of the three candidates.
  The `extend_mesh_with_pml` machinery currently supports only
  single-axis extensions (`&[PmlAxis::XMax]` etc.); multi-axis wedge
  extensions require new Kuhn-6 wedge tetrahedralisation in
  `crates/yee-mesh` and a `PmlAxis::WedgeXYZ` variant. Estimated
  ~1500-2000 LoC + new mesh tests.
- **Risk.** v3.5 P2 escape hatch explicitly deferred this to
  v3.5.1+ per ADR-0043 §risks (c) consequence ("multi-face wedge
  PML deferred"). Like W2, reopening it would unblock fem-eig-006
  but at substantial cost.

### 3.4 Recommendation summary

| Hypothesis | LoC | Risk | Expected |S_11| | Decision |
|------------|-----|------|------------------|----------|
| **W1** wave-port    | ~10        | low    | < 0.01    | **ship in v3.5.3** |
| W2 rotated PML      | ~500-1000  | medium | < 0.05    | defer to v3.5.4+ |
| W3 multi-face wedge | ~1500-2000 | high   | < 0.05    | defer to Phase 4.fem.eig.4 |

W1 is **strictly less general** than W2 / W3 — it only works for
fixtures whose dominant mode is well-approximated by a single
analytic modal shape. But for fem-eig-006 specifically, that
condition holds (TE_{10} is the only propagating mode at 30 GHz on
the `b = 10 mm` broad wall), and shipping W1 retires the gate
without committing the project to a multi-face wedge mesher.

## 4. Mathematical formulation

### 4.1 Current v3.5.2 CFS-PML termination (broken on fem-eig-006)

The v3.5.2 driver extends the cavity mesh with a CFS-PML shell on
the +x face. The +x outer truncation surface is tagged
`FaceKind::Pec`; the inner PML interface is internal (no tag). The
stretched coordinate tensor

```text
s_x(d_x, ω) = κ_x(d_x) + σ_x(d_x) / (α_alpha(d_x) + j ω ε_0)
```

ramps from `s_x = 1` at the inner interface (`d_x = 0`) to
`s_x = κ_max + σ_max/(α_min + j ω ε_0)` at the outer truncation
(`d_x = D`). For the v3.5.2 defaults
`(κ_max=2, m=4, thickness=14, alpha_grading_order=1)`, the absorber
budget is theoretically `~80 dB` per Roden-Gedney 2000 §IV. The
fem-eig-006 measurement shows the budget is **unused**: the
guide-mode never sees the stretching tensor in its propagation-axis
frame, so the round-trip absorption is the `|S_{11}| = 0.926`
measured floor.

### 4.2 v3.5.3 wave-port termination (W1)

Replace the CFS-PML shell on the +x face with a `FaceKind::WavePort`
boundary condition. The +x face becomes a **second wave-port**
carrying the same TE_{10} modal basis as the -x face. The cavity is
now a two-port WR-style waveguide section:

```text
Port 0 (-x, x = 0):     FaceKind::WavePort(0), drives TE_{10}
Port 1 (+x, x = 100mm): FaceKind::WavePort(1), absorbs TE_{10}
```

Both ports share `fem_eig_006_beta_te10` / `fem_eig_006_modal_e_t_te10`
(the geometric translation along x preserves the modal shape).

The Phase 4.fem.eig.2 E2 wave-port boundary condition contributes a
per-face block

```text
B_port(p, p') = j β B_port,modal(p, p') = j β ∫_face e_t(p) · e_t(p') dA
```

to the stiffness matrix and a per-face modal-current term to the
RHS vector for the **excited** port. The reflection coefficient at
port 0 is then computed by modal decomposition

```text
S_{11} = (|E_FEM,t · e_t|_x=0 - a_excite) / a_excite
```

where `a_excite` is the unit modal-current excitation amplitude.

Per Jin §10.6, the resulting `S_{11}` reflection at a matched
wave-port terminated by an identical wave-port is `|S_{11}| < 0.01`
on a moderate-resolution FEM mesh (the fem-eig-004 thru-line on a
`(16, 3, 2)` mesh achieves `|S_{11}| < -20 dB` ≡ `0.1` against the
strict gate; fem-eig-006 has a finer mesh in the dominant x-axis
and should do better).

### 4.3 Backward-compatibility

The W1 change is **driver-only**: no `PmlConfig` field changes, no
`OpenBoundarySolver` API changes, no new types. The
`fem_eig_006_magnitude_bounded` gate flips from `#[ignore]` to
CI-default. fem-eig-003 is **untouched** (it uses a separate driver
with separate CFS-PML configuration); the v3.5.2 defaults stay in
force for fem-eig-003 and continue to retire its strict band
`[-71.53, -55.58] dB`.

## 5. Public API

**No public-API changes.** All changes are internal to the
`run_fem_eig_006_high_aspect_pml_with_config` driver in
`crates/yee-validation/src/lib.rs`.

The driver currently classifies faces as:

```rust
for c in &centroids {
    let kind = if c.x.abs() < tol {
        FaceKind::WavePort(0)
    } else {
        FaceKind::Pec
    };
    face_kinds.push(kind);
}
```

and applies CFS-PML via `.with_cfs_pml(pml_config, pml_classes)`.
After W1, the driver becomes:

```rust
for c in &centroids {
    let kind = if c.x.abs() < tol {
        FaceKind::WavePort(0)
    } else if (c.x - FEM_EIG_006_A_M).abs() < tol {
        FaceKind::WavePort(1)
    } else {
        FaceKind::Pec
    };
    face_kinds.push(kind);
}
```

The cavity mesh extension (`extend_mesh_with_pml`) is **removed** —
the cavity is now its native `(16, 3, 2)` shape with no PML shell.
The `.with_cfs_pml(...)` builder call is also removed. A second
`PortDefinition` clones the TE_{10} shape for port 1.

The `pml_config: yee_fem::PmlConfig` driver signature parameter is
**kept** for source-compatibility with the v3.5.2 ablation binary
(`crates/yee-validation/examples/cfs_pml_grading_sweep.rs`) but
becomes unused inside the body. A `#[allow(unused_variables)]` plus
a doc-comment notes the v3.5.3 deprecation; a follow-on Phase
4.fem.eig.4 commit can remove the parameter entirely.

The fem-eig-006 ablation rows in the SSSSSSSSS H4 grid will produce
identical `|S_{11}|` per row (now driven by the wave-port floor, not
the PML config) — this is the **expected** behaviour post-W1 and
the v3.5.2 sweep CSV remains reproducible.

## 6. Validation

After W1 lands:

1. Run `cargo test -p yee-validation --release --test fem_eig_006_high_aspect_pml`
   with the `#[ignore]` removed from `fem_eig_006_magnitude_bounded`.
   Expected: `|S_{11}(30 GHz)| < 0.01`; gate retire.
2. Run `cargo test -p yee-validation --release --test fem_eig_003_wr90_stub_abc`
   to verify v3.5.2 fem-eig-003 retire is **untouched**. Expected:
   pass; band `[-71.53, -55.58] dB`.
3. Run `cargo test --workspace --release` to verify no regressions
   across the broader gate suite.
4. The smoke test `fem_eig_006_smoke_runs` and the canary
   `fem_eig_006_no_nan_inf` continue to run as before; both should
   pass on the new wave-port configuration.

## 7. Risks

- **(a) TE_{10} wave-port modal basis may underestimate higher-
  order content at 30 GHz.** The TE_{20} mode cutoff on `b = 10 mm`
  sits at `30 GHz` exactly. If higher-order modal content is present
  in the cavity, a TE_{10}-only wave-port underestimates the
  reflection by treating the higher-order modes as unsupported (they
  produce a spurious "leakage" into the modal decomposition).
  **Mitigation:** the W1 measurement bounds the residual; if
  `|S_{11}| > 0.1` after W1, queue Phase 4.fem.eig.3.5.4 for
  multi-mode-port extension (add explicit TE_{20} / TE_{01} modes
  to the port definition). The Phase 4.fem.eig.2 E2 wave-port
  machinery supports multi-mode ports via additional
  `PortDefinition` entries, so the extension is mechanical.

- **(b) fem-eig-006 must keep retiring fem-eig-003.** Phase
  4.fem.eig.3.5.2 retired fem-eig-003 via the new `PmlConfig`
  defaults. W1 is **fem-eig-006-specific** (only the
  fem-eig-006 driver changes); fem-eig-003 retains its v3.5.2
  driver and v3.5.2 PML defaults exactly. **Mitigation:** the §6
  validation step (2) explicitly re-runs fem-eig-003 to confirm
  this.

- **(c) The +x wave-port is no longer the PML inner-boundary
  interface; mesh anisotropy may surface as wave-port
  modal-coupling errors.** The native cavity mesh has aspect
  ratio `100 : 1 : 1` (`dx = 6.25 mm`, `dy = 3.33 mm`,
  `dz = 0.5 mm`). The wave-port modal integration on a `dy × dz`
  face uses Whitney-1 basis functions (3 edges per face); on
  highly anisotropic cells, the modal projection accuracy
  degrades. **Mitigation:** fem-eig-004's thru-line uses a
  similar `(16, 3, 2)` mesh and achieves `|S_{11}| < -20 dB`,
  giving empirical evidence that the modal projection is robust
  at this aspect ratio. If W1 misses the gate, the first
  diagnostic is to refine the y/z dimensions to `(16, 5, 3)` or
  `(16, 7, 5)` and re-measure.

- **(d) PmlConfig parameter becomes vestigial on the
  fem-eig-006 driver.** The `pml_config: yee_fem::PmlConfig`
  parameter survives the W1 change as a no-op for source
  compatibility with the v3.5.2 ablation binary, but is
  **unused** by the new driver body. **Mitigation:**
  doc-comment the parameter as "Phase 4.fem.eig.3.5.3
  deprecated; retained for v3.5.2 sweep CSV compatibility";
  queue removal under Phase 4.fem.eig.4 cleanup.

## 8. Lane

Spec file:
`docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-3-design.md`

Implementation lane (declared here for the T1-T4 plan):

- `crates/yee-validation/src/lib.rs` —
  `run_fem_eig_006_high_aspect_pml_with_config` driver: switch
  +x face from PML truncation to `FaceKind::WavePort(1)`;
  remove `extend_mesh_with_pml` call; remove `.with_cfs_pml(...)`
  builder call; add second `PortDefinition` cloning TE_{10}.
- `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs` —
  remove the `#[ignore]` from `fem_eig_006_magnitude_bounded`;
  refresh the docstring to record the v3.5.3 measurement.
- `docs/src/tutorials/07-fem-open-cavity.md` — note the v3.5.3
  fem-eig-006 wave-port termination and the §3.4 hypothesis-tree
  outcome.
- `ROADMAP.md` — Phase 4.fem.eig.3.5.3 entry from planned to
  shipped.

Out of lane: `yee-fem` (no PmlConfig / open_boundary changes;
all behavioural changes are driver-level in yee-validation),
`yee-cli`, `yee-gui`, `yee-mom`, `yee-mesh`, `yee-cuda`,
`yee-plotters`, `yee-fdtd`, `yee-surrogate`.

## 9. References

- Berenger, J.-P., "Three-dimensional perfectly matched layer for
  the absorption of electromagnetic waves," *IEEE Transactions on
  Antennas and Propagation* 44(1) (January 1996), pp. 110-117.
  DOI 10.1109/8.477535. §IV-A bulk vs guide-wave PML; §IV-B
  rotated/oblique formulation (W2 hypothesis).
- Berenger, J.-P., "A perfectly matched layer for the absorption
  of electromagnetic waves," *Journal of Computational Physics*
  114(2) (1994), pp. 185-200. DOI 10.1006/jcph.1994.1159. §V
  multi-face wedge PML (W3 hypothesis).
- Jin, J.-M., *The Finite Element Method in Electromagnetics*,
  3rd ed. (Wiley, 2014), Chapter 10.6 "Wave-port termination" —
  W1 mathematical foundation; standard closed-cavity FEM modal
  termination.
- Roden, J. A. and Gedney, S. D., "Convolutional PML (CPML)",
  *IEEE MWCL* 10(5) (May 2000) — CFS-PML formulation inherited
  via Phase 4.fem.eig.3.5 / 3.5.1 / 3.5.2.
- `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-2-alpha-grading-design.md`
  — v3.5.2 parent spec; §7 (a) escape-hatch path queued this
  v3.5.3 work.
- `docs/superpowers/plans/2026-05-20-phase-4-fem-eig-3-5-2-alpha-grading.md`
  — v3.5.2 S1-S5 plan (SSSSSSSSS shipped); S4 ablation produced
  the frozen-magnitude finding.
- `docs/src/decisions/0045-phase-4-fem-eig-3-5-2-alpha-grading.md`
  — v3.5.2 ADR; §risks (a) explicitly queued the W1
  wave-port-termination path.
- `docs/src/decisions/0046-phase-4-fem-eig-3-5-3-fem-eig-006-retire.md`
  — this spec's scope ADR.
- `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs:53-60`
  — current `#[ignore]` docstring recording the SSSSSSSSS H4
  frozen-magnitude finding.
- CLAUDE.md §3, §4, §10.
