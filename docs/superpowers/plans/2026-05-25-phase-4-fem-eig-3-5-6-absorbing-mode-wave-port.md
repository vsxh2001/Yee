# Phase 4.fem.eig.3.5.6 — Lee-Mittra Absorbing-Mode Wave-Port — Implementation Plan

**Companion spec:** `docs/superpowers/specs/2026-05-25-phase-4-fem-eig-3-5-6-absorbing-mode-wave-port-design.md`  
**ADR:** ADR-0070  
**Lane:** `crates/yee-fem/src/**`, `crates/yee-validation/src/lib.rs`,
          `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`  
**PATTERN FILE:** `crates/yee-fem/src/element.rs::assemble_port_face_block_gauss_pts`
  (imitate this function for the new projected variant)  
**VERIFICATION COMMAND:**
```
cargo test -p yee-validation --test fem_eig_006_high_aspect_pml --release -- --include-ignored
```
Expected: all 3 tests green (including `fem_eig_006_magnitude_bounded` if gate passes).
Fallback (escape hatch): 2 green (smoke + finite), 1 ignored with updated notes.

**ESCAPE HATCH:** blocked > 15 min → surface and stop.

---

## Step L1 — Element-layer: `assemble_port_face_block_projected_gauss_pts`

**File:** `crates/yee-fem/src/element.rs`  
**After:** `assemble_port_face_block_gauss_pts` (around line 890)

Add a new `pub fn assemble_port_face_block_projected_gauss_pts`:

```rust
/// Rank-1 modal-projected wave-port face block (Lee-Mittra 1997 §IV).
///
/// Computes
///
/// ```text
/// block[i,j] = j · β_eff / μ_r · Σ_g w_g · [(n̂ × N_i(ξ_g)) · e_t_g]
///                                            · [(e_t_g · n̂ × N_j(ξ_g))]
/// ```
///
/// using the exact Whitney-1 basis at the same 3-point Gauss set as
/// [`assemble_port_face_block_gauss_pts`]. The result is a **rank-1**
/// matrix in the edge-DoF space (outer product of the modal projections).
///
/// When `β_eff = 0` the block vanishes identically (same as the scalar
/// path for evanescent modes — no special-case needed).
pub fn assemble_port_face_block_projected_gauss_pts(
    face_vertices: [Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
    beta_eff: f64,
    modal_e_t_at_gauss_pts: [Vector3<f64>; 3],
    mu_r_face: f64,
) -> SMatrix<Complex64, 3, 3>
```

Implementation follows `assemble_port_face_block_gauss_pts` step by step, but the
inner loop accumulates `w_g * a_i * a_j` instead of `w_g * (n̂ × N_i) · (n̂ × N_j)`,
where `a_i = (n̂ × N_i(ξ_g)) · e_t_g`.

Unit tests (in a `#[cfg(test)] mod tests_projected_block` below the function):

- `L1_rank1_structure`: for a non-zero e_t, verify `block[0,1] * block[2,0] ≈ block[0,0] * block[2,1]`
  (rank-1 outer-product identity `B[i,j] × B[k,l] = B[i,l] × B[k,j]`).
- `L1_orthogonal_e_t_zero`: when `e_t_g = (0,0,1)` on a face in the xy-plane (n̂ = (0,0,1)),
  n̂ × N is in-plane so the dot with out-of-plane e_t is zero → block is zero.
- `L1_uniform_e_t_matches_scalar`: when `e_t_g` is aligned with `(n̂ × N_avg_normed)` for
  all Gauss points, the projected block should match the scalar block up to the projection
  factor (smoke check only, can use an analytic simple triangle).

Estimated LOC: ~120 (60 implementation + 60 tests).  
Verification: `cargo test -p yee-fem --lib -- element::tests_projected_block` exits 0 in < 5 s.

---

## Step L2 — Element-layer: `assemble_port_face_block_projected` (centroid path)

**File:** `crates/yee-fem/src/element.rs`  
**After:** `assemble_port_face_block` (around line 598)

Add the centroid-approximation variant mirroring `assemble_port_face_block`:

```rust
/// Centroid-approximation variant of
/// [`assemble_port_face_block_projected_gauss_pts`] for the
/// non-coupled-Whitney path.
///
/// ```text
/// block[i,j] = j · β_eff / μ_r · face_area · [(n̂ × t_i) · e_t_c]
///                                             · [(e_t_c · n̂ × t_j)]
/// ```
pub fn assemble_port_face_block_projected(
    face_vertices: [Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
    beta_eff: f64,
    modal_e_t_at_centroid: Vector3<f64>,
    mu_r_face: f64,
) -> SMatrix<Complex64, 3, 3>
```

Unit tests parallel to L1 (same structure / orthogonal / smoke checks). Reuse the analytic
triangle from L1.

Estimated LOC: ~80 (40 + 40 tests).  
Verification: `cargo test -p yee-fem --lib -- element::tests_projected_centroid` exits 0.

---

## Step L3 — `PortDefinition`: add `absorbing_complement` field + builder

**File:** `crates/yee-fem/src/open_boundary.rs`  
**Section:** `PortDefinition` struct (around line 616)

```rust
pub struct PortDefinition {
    pub modes: Vec<PortMode>,
    /// Enable Lee-Mittra first-order absorbing-mode complement
    /// (Phase 4.fem.eig.3.5.6). Default `false` → existing scalar
    /// wave-port stiffness (backward-compat).
    pub absorbing_complement: bool,
}
```

Add builder:
```rust
impl PortDefinition {
    /// Enable the Lee-Mittra first-order absorbing-mode complement
    /// on this port face. See spec §3.2.
    pub fn with_absorbing_complement(mut self) -> Self {
        self.absorbing_complement = true;
        self
    }
}
```

Ensure `PortDefinition::single_mode` (and any other constructors) leaves `absorbing_complement: false`.

Unit test: `port_definition_default_absorbing_complement_is_false` — check `PortDefinition::single_mode(...).absorbing_complement == false`.

Estimated LOC: ~20.  
Verification: `cargo test -p yee-fem --lib -- open_boundary` exits 0.

---

## Step L4 — `scatter_port_face_gauss`: Lee-Mittra branch

**File:** `crates/yee-fem/src/open_boundary.rs`  
**Function:** `scatter_port_face_gauss` (line ~2173)

Add imports at the top of the function:
```rust
use crate::element::assemble_port_face_block_projected_gauss_pts;
const C0: f64 = 2.997_924_58e8;
```
(C0 may already be in scope — check the existing import.)

After the function signature, add an early branch on `port.absorbing_complement`:

```rust
if port.absorbing_complement {
    // Lee-Mittra first-order absorbing-mode BC (spec §3.2):
    // K = jk₀ B_face + ∑_m j(β_m − k₀) R_m
    // = assemble_port_face_block_gauss_pts(..., k0) + ∑_m proj_block(β_m − k₀, e_t_m)
    let k0 = omega / C0;
    let k_full = assemble_port_face_block_gauss_pts(
        face_vertices, face.normal, Complex64::new(k0, 0.0), 1.0
    );
    let mut k_lee_mittra = k_full;
    for mode in &port.modes {
        let beta = (mode.beta_mode)(omega);
        let beta_eff = beta - k0;
        let e_t_gauss = /* sample modal_e_t at the 3 Gauss points (same as existing) */;
        let r_m = assemble_port_face_block_projected_gauss_pts(
            face_vertices, face.normal, beta_eff, e_t_gauss, 1.0
        );
        k_lee_mittra += r_m;
    }
    // Scatter k_lee_mittra into triplets (same PEC-precedence guard as existing).
    // Then scatter RHS (unchanged — same a_inc × 2jβ_m loop).
    // ... (mirror existing scatter loop)
    return;
}
// Existing non-absorbing_complement path below (unchanged)
for mode in &port.modes { ... }
```

Keep the existing loop BELOW the branch — it runs only when `absorbing_complement = false`.

**Critical:** RHS accumulation (lines after stiffness scatter in existing loop) must be
preserved verbatim in the absorbing_complement branch's RHS sub-loop.

Estimated LOC: ~60 (new branch + reuse existing RHS loop).  
Verification: `cargo test -p yee-fem --lib -- open_boundary::tests` exits 0 in < 30 s.

---

## Step L5 — `scatter_port_face`: Lee-Mittra branch (centroid path)

**File:** `crates/yee-fem/src/open_boundary.rs`  
**Function:** `scatter_port_face` (line ~2099)

Mirror L4 using `assemble_port_face_block_projected` (centroid approximation) and
`assemble_abc_face_block(face_vertices, face.normal, k0, 1.0)` for the k₀ full block.

Estimated LOC: ~50.  
Verification: same as L4.

---

## Step V1 — Validation fixture: enable `absorbing_complement` on port_1

**File:** `crates/yee-validation/src/lib.rs`  
**Function:** `run_fem_eig_006_high_aspect_pml_with_config` (line ~4030)

Change port_1 construction to:
```rust
let port_1 = PortDefinition {
    modes: vec![te10_mode, te20_mode, te01_mode],
    absorbing_complement: true,
}.with_absorbing_complement(); // or just set the field directly
```

No other change to the driver.

Estimated LOC: ~3.  
Verification: `cargo test -p yee-validation --test fem_eig_006_high_aspect_pml --release fem_eig_006_no_nan_inf` exits 0 first (smoke check), then run magnitude gate.

---

## Step V2 — Gate disposition

After measuring `|S₁₁|` with the absorbing complement:

**If `|S₁₁| < 0.1`:**
1. Un-ignore `fem_eig_006_magnitude_bounded` in
   `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`.
2. Update the test docstring with the new measurement.
3. Update `run_fem_eig_006()` in `lib.rs` to report `CaseStatus::Passed`.

**If `|S₁₁| ≥ 0.1`:**
1. Keep `#[ignore]` on `fem_eig_006_magnitude_bounded`.
2. Update the ignore reason with: new measurement, % improvement from 0.955, and escape-hatch
   diagnosis (e.g., "Lee-Mittra first-order complement improves |S₁₁| from 0.955 to X.XX;
   higher-order absorbing BC queued for Phase 4.fem.eig.3.5.7").
3. Surface the finding in the agent report.

---

## Step D1 — ADR, ROADMAP, CLAUDE.md update

**Files:** `docs/src/decisions/0070-...md`, `docs/src/SUMMARY.md`, `ROADMAP.md`, `CLAUDE.md`

- ADR-0070 records the decision, measurements, and outcome.
- ROADMAP status snapshot: "Phase 4.fem.eig.3.5.6 shipped" (or "attempted; finding: X").
- CLAUDE.md §10: update fem-eig-006 status line.

---

## Verification sequence

```bash
# L1+L2: element unit tests
cargo test -p yee-fem --lib -- tests_projected 2>&1 | tail -5

# L3: PortDefinition unit test
cargo test -p yee-fem --lib -- open_boundary 2>&1 | tail -5

# Regression guards (must not regress)
cargo test -p yee-validation --release fem_eig_001 fem_eig_002 fem_eig_004 fem_eig_005 2>&1 | tail -10

# Gate measurement
cargo test -p yee-validation --test fem_eig_006_high_aspect_pml --release -- --include-ignored 2>&1 | tail -10

# Lint (required before reporting done)
cargo clippy -p yee-fem -p yee-validation --all-targets -- -D warnings 2>&1 | tail -5
cargo fmt --check -p yee-fem -p yee-validation 2>&1 | tail -5
```
