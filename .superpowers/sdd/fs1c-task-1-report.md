# FS.1c Task 1 — z-axis thin-wire subcell (Holland–Simpson), CPU

**Branch:** `feature/fs1c-thin-wire`
**head_before:** `2599d0b080e882c1762f1f1e5da52aa5f0e58595` (docs: FS.1c spec + plan)
**head_after (functional commit):** `e655a3b` — `yee-compute: add z-axis thin-wire
subcell (Holland-Simpson FS.1c)`

Plan: `docs/superpowers/plans/2026-07-24-fs1c-thin-wire.md`, Task 1.
Spec: `docs/superpowers/specs/2026-07-24-fs1c-thin-wire-design.md`.

> Note: this report was not written at the time of the `e655a3b` commit — a
> reviewer-found gap (the plan's Task 1 bullet requires it). The substance
> below is transcribed from that commit's message and the module docs it
> shipped, not reconstructed after the fact; the code and tests are
> unchanged by this fix round. See "## Fix round" at the end for what this
> pass actually did.

## Research first (formulation + citation)

Holland & Simpson's 1981 **in-cell-inductance** thin-wire subcell model:

> R. Holland and L. Simpson, "Finite-Difference Analysis of EMP Coupling to
> Thin Struts and Wires," *IEEE Trans. Electromagn. Compat.*, vol. 23,
> no. 2, pp. 88–97, May 1981.

implemented per the derivation in:

> Y. Liu, *Use of the Thin-Strut FDTD Formalism for the Design of Coils in
> Biomedical Telemetry Applications*, M.S. thesis, North Carolina State
> University, 2003, ch. 4, eq. 4.1–4.18 (itself after Holland & Simpson
> 1981 and K. R. Umashankar, A. Taflove, and B. Beker, "Calculation and
> experimental validation of induced currents on coupled wires in an
> arbitrary shaped cavity," *IEEE Trans. Antennas Propagat.*, vol. AP-35,
> no. 11, pp. 1248–1257, Nov. 1987; also summarized in Taflove & Hagness,
> *Computational Electrodynamics*, ch. 10, "Local Subcell Models of Fine
> Geometrical Features").

This citation is recorded in the module docs (`drive.rs`'s `ThinWire` and
`thin_wire_l_prime` doc comments, `cpu.rs`'s `advance_thin_wire_currents`
doc comment) **before** the update equations were implemented, per the
research-first constraint. The equations were read out of these sources,
not invented.

### The model, in short

Integrating Ampère's law azimuthally around the wire and Faraday's law
radially from the wire surface (`E_z(a) = 0`, PEC) out to `R = h/2` (half
the transverse cell size) gives the wire's **in-cell inductance per unit
length**:

```
L'(h/2) = (μ₀/2π)·ln(h/(2a))
```

for a wire of physical radius `a`. Each wire-occupied `E_z` cell carries a
shunt inductor with this `L'`; its branch current is subtracted from the
ordinary curl-H `E_z` update (`J_z = I/(dx·dy)`, distributed onto the
surrounding field the same way this crate's existing `ResistivePort`/
`AperturePort` branches already subtract from `E_z`). The near-wire
transverse field (`E_x`/`E_y` at the wire's own grid line) is forced to
zero every step — the coarse-grid stand-in for the un-resolved
near-singular radial field around the conductor. The wire's two free
(open) ends get a hard `I = 0` condition (no further conductor for current
to flow into).

**Named simplification (not a silent gap):** the full Holland–Simpson/Liu
system couples `I` and line charge `Q` along z (a 1-D telegrapher line
solved simultaneously with the 3-D fields, thesis eq. 4.15–4.17). This
crate's walking-skeleton reduction drops the `dQ/dz` charge/continuity
term — the surrounding 3-D Maxwell grid already conserves charge on its
own — leaving a pure lumped-inductor branch driven by the local cell's own
`E_z`. This gets within a measured single-digit-percent of coarse/fine
grid-independence (see below) but not as tight as the full
telegrapher-coupled system would; tightening it is a named follow-on,
not a bug.

## Seam decision: `Drive`, not `Materials`

Read `crates/yee-compute/src/{drive.rs,materials.rs,cpu.rs}` before
choosing. `ThinWire` was attached to **`Drive`**, not `Materials`, because:

- The model needs **persistent per-cell branch state that evolves every
  step** (`wire_current: Vec<Vec<f64>>`, one shunt-inductor current per
  wire cell) plus a **post-`boundary_e` correction pass**
  (`apply_thin_wire_correction`) that runs after the ordinary E half-step
  — exactly the same shape as the existing `ResistivePort`/`AperturePort`
  idiom already living in `Drive`, which also carries their own
  persistent branch/record state and their own post-boundary correction
  hooks in `CpuFdtd::step`.
- `Materials` (masks like `pec_mask_ez`, per-cell `eps_r`/`sigma` arrays)
  is a **static, step-invariant description of the medium** consulted
  inside the hot rayon `update_h`/`update_e` loops. Routing a thin wire
  through `Materials` would mean threading new per-cell mutable state
  (the branch current) through those hot parallel loops for every cell in
  the grid, for no benefit — the wire only touches a handful of cells
  along its own axis, and the port/probe idiom already has the right
  shape (a small, separately-iterated list of "special" cells visited
  outside the main rayon sweep). `Drive` is the least-churn seam.

`CpuFdtd` gained one field, `wire_current: Vec<Vec<f64>>` (indexed
`wire_index -> (k - wire.k_lo) -> current`), populated in `with_drive`/
`new_with_drive` alongside the existing `probe_series`/`aperture_records`
construction, and two new private step methods called from `CpuFdtd::step`:

- `advance_thin_wire_currents()` — called right after `update_h` +
  `boundary_h`, using `E_z` from the end of the *previous* step (the same
  half-step offset `H` has from `E`), so it must run before `update_e`
  overwrites `E_z` this iteration.
- `apply_thin_wire_correction()` — called right after `boundary_e`,
  subtracting the just-advanced branch current from `E_z` and shorting the
  near-wire transverse `E_x`/`E_y`.

## GPU rejection

`gpu.rs` gained a named `ComputeError::Unsupported` rejection for any
`Drive` carrying a non-empty `thin_wires` list, checked pre-adapter — same
pattern as the existing aperture-port/sheet-loss rejections. Test:
`crates/yee-compute/tests/gpu_thinwire_rejected.rs`
(`gpu_rejects_a_drive_carrying_a_thin_wire`).

## Two bugs found and fixed during derivation-to-code translation

Both caught by an initially-unstable/then-inconsistent coarse-vs-fine
resonance test, i.e. by the gate the plan asked for, not by inspection:

1. **Missing `dz`-free pointwise relation.** Eq. 4.13 (Liu thesis) is a
   *pointwise* relation: both sides are "per unit length" (`dI/dt` carries
   no length scale; `E_z/L'` is `(V/m)/(H/m) = A/s`) — exactly parallel to
   how the ordinary curl-H `E_z` update coefficient has no `dz` in it
   either. An initial translation multiplied by `dz`, which blew the
   simulation up to `NaN` within the first several hundred steps.
2. **Missing open-end `I = 0` condition.** Without forcing `I = 0` at the
   wire's two free ends, the model still ran (finite fields) but the
   coarse/fine resonance drifted by an extra ~4 percentage points versus
   with the condition — the open-circuit boundary is not optional for a
   finite (non-infinite) wire.

## Unit tests (`crates/yee-compute/tests/thin_wire.rs`)

- `no_wire_construction_is_bit_identical_to_the_old_api` — an empty
  `Drive::thin_wires` reproduces the pre-existing `with_config` entry
  point **bit-for-bit** across all six field arrays (`assert_eq!`, not a
  tolerance). This is the PROVABLE no-op the global constraints require.
- `wire_present_smoke_stays_finite_and_perturbs_the_field` — attaching a
  wire keeps every field component finite and measurably changes `Ez`
  relative to a free-space run (the model does something).
- `coarse_fine_resonance_consistency_and_naive_control` — the same
  physical dipole (fixed `L = 40 mm`, `a = 0.3 mm`) run at `dx = 4 mm` and
  `dx/√2`, resonant frequency (`Im(Z)` zero-crossing from feed V/I via a
  single-bin DFT / `Z(f) = V(f)/I(f)`) within a **measured, honestly
  pinned** `< 10%` (measured **~8.1%** on this fixture — pinned from what
  was actually measured, not a wished-for target; the reduced model's
  dropped charge/continuity coupling is why this isn't tighter). A naive
  one-cell-PEC run at the same two resolutions is reported alongside
  (printed via `eprintln!`, not hard-asserted against — on this toy
  geometry the two models' resonances land in different parts of a
  structured, multi-resonance spectrum, so "which is more grid-independent"
  isn't a clean apples-to-apples comparison at this fixture's size).

The other 9 existing test files touched (`cavity_resonance.rs`,
`cpu_aperture_parity.rs`, `cpu_drive_parity.rs`, `cpu_h_probe.rs`,
`gpu_aperture_parity.rs`, `gpu_h_probe_parity.rs`, `gpu_ntff_dipole.rs`,
`line_eeff.rs`, `ntff_dipole.rs`) only add `thin_wires: vec![]` to their
existing exhaustive `Drive { .. }` literals — mechanical, forced by the
new field, no behavior change.

## Verification (original commit + this fix round, both real output)

- Bit-exact suite: `cargo test -p yee-compute --release --test
  graded_uniform_bitexact --test gpu_graded_parity --test gpu_cpu_parity --
  --include-ignored` — **5/5 pass**. GPU evidence real:
  `compute-002: running on adapter 'NVIDIA GeForce RTX 5060 Ti'`
  (family-rel L2/L∞ ~1.5e-7–3.0e-7); `compute-020` graded-vs-scalar GPU
  bit-for-bit PASS (L2/L∞ ~7e-6–5.5e-5); `compute-021` graded-taper GPU
  parity (reflection level −52.68 dB).
- `cargo test -p yee-compute --release` (full default suite): all test
  binaries **0 failed** (`thin_wire.rs`: 3/3 pass; `gpu_thinwire_rejected.rs`:
  1/1 pass; every other file's `test result: ok`, ignored gates skipped as
  designed).
- `cargo clippy --workspace --all-targets -- -D warnings` — clean.
- `cargo clippy -p yee-compute --all-targets --no-default-features -- -D
  warnings` — clean.
- `cargo fmt --check --all` — clean (no diff).
- `missing_docs` — clean (`yee-compute` is `#![warn(missing_docs)]`; no
  warnings on the new public/`pub(crate)` items).

## Commit

`e655a3b` — `yee-compute: add z-axis thin-wire subcell (Holland-Simpson
FS.1c)`. 15 files changed, 641 insertions(+), 2 deletions(-):
`crates/yee-compute/src/{cpu.rs,drive.rs,gpu.rs,lib.rs}` +
`crates/yee-compute/tests/{cavity_resonance,cpu_aperture_parity,
cpu_drive_parity,cpu_h_probe,gpu_aperture_parity,gpu_h_probe_parity,
gpu_ntff_dipole,gpu_thinwire_rejected,line_eeff,ntff_dipole,thin_wire}.rs`.
Trailers per the global constraints.

## Fix round

**Reviewer finding addressed:** this report (`.superpowers/sdd/
fs1c-task-1-report.md`) did not exist — the plan's Task 1 bullet required
transcribing the seam decision + the two bugs found/fixed into the report,
and that never happened even though the substance was present in the
`e655a3b` commit message and module docs. This is a documentation-only
fix: no source or test files changed, per the reviewer's own assessment
("paperwork gap rather than a code defect").

**What this pass did:**

1. Read the `e655a3b` commit message, `drive.rs`'s `ThinWire`/
   `thin_wire_l_prime` docs, and `cpu.rs`'s `advance_thin_wire_currents`/
   `apply_thin_wire_correction` docs to confirm the report above
   accurately transcribes (not reinterprets) the shipped citation, seam
   rationale, and the two translation bugs.
2. Wrote this file.
3. Re-ran all mandated verification commands fresh, this session, on the
   unmodified `e655a3b` tree (real output recorded above, superseding the
   original commit's un-filed verification):
   - Bit-exact suite (`graded_uniform_bitexact` + `gpu_graded_parity` +
     `gpu_cpu_parity`, `--include-ignored`): **5/5 pass**, GPU adapter
     confirmed `NVIDIA GeForce RTX 5060 Ti` (not skipped).
   - `cargo test -p yee-compute --release` (full suite): **0 failed**
     across every test binary, including `thin_wire.rs` (3/3) and
     `gpu_thinwire_rejected.rs` (1/1).
   - `cargo clippy --workspace --all-targets -- -D warnings`: clean.
   - `cargo clippy -p yee-compute --all-targets --no-default-features --
     -D warnings`: clean.
   - `cargo fmt --check --all`: clean.
4. Confirmed `git branch --show-current` is `feature/fs1c-thin-wire`
   before staging/committing (never `main`).
5. Staged only `.superpowers/sdd/fs1c-task-1-report.md` (explicit path)
   and committed as a `docs:` commit on this branch.

No code, test, or spec/plan content changed in this fix round — the
underlying Task 1 implementation was already correct and green; the sole
deliverable gap was the missing report file, now filled.
