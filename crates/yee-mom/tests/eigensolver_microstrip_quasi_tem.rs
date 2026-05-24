//! Phase 1.3.1.2 validation gate — quasi-TEM microstrip cross-section
//! eigensolve vs the Hammerstad-Jensen `ε_eff`.
//!
//! The cross-section eigensolver's purpose is **quasi-TEM microstrip
//! wave-ports**, but its closed-guide dominant-mode selection
//! (`solve_dense_mixed`) gathers candidates from the *cutoff* pencil and
//! floors out the `k_c² ≈ 0` cluster, discarding the quasi-TEM mode (which
//! also sits at `k_c² ≈ 0`) along with the curl-free gradient nulls — so on
//! a microstrip it returns the box's air-region mode (`ε_eff ≈ 1`), not the
//! quasi-TEM (ADR-0059). Phase 1.3.1.2 adds a **separate quasi-TEM selection
//! path** (`solve_dense_mixed_quasi_tem`, reached here through the
//! `__internal::quasi_tem_beta_for_test` test surface) that targets the
//! quasi-TEM by a TEM-scale β-direct shift-invert ladder and discriminates
//! it from the gradient nulls with the converged-eigenvector transverse
//! screen.
//!
//! ## What this gate checks (DoD-2)
//!
//! On a **canonical shielded 50Ω FR-4 microstrip** (`ε_r = 4.4`, `w/h = 1.9`)
//! the quasi-TEM `ε_eff = (β/k₀)²` reconciles to within a **loose** tolerance
//! of the **Hammerstad-Jensen open-line `ε_eff`** (Balanis *Antenna Theory*
//! 3rd ed. Eq. 14-1; Pozar §3.8) — the published microstrip benchmark. The
//! tolerance is loose because the **shielding box truncates the open
//! geometry**: a finite box forces some field that an open line carries far
//! above the strip back into the substrate/air region, perturbing `ε_eff`.
//! A larger box tightens the agreement; the box here is sized so the shift
//! is within the loose band (see [`BOX_W_OVER_STRIP`] / [`BOX_H_OVER_SUB`]).
//!
//! ## How the microstrip is modeled (shielded, strip-as-hole)
//!
//! A microstrip is a **two-conductor** line (signal strip + ground/shield);
//! a quasi-TEM mode exists only because there are two separated conductors
//! (a single hollow conductor supports no TEM mode). The cross-section
//! eigensolver imposes PEC Dirichlet on every **mesh-boundary** edge (an
//! edge in exactly one triangle). We exploit that: the signal strip is a
//! rectangular **hole** in the mesh (its cells are not triangulated), so the
//! hole's border edges are mesh-boundary = PEC — a second PEC conductor
//! inside the outer PEC box (the shield + ground plane). Substrate
//! (`y < h`, `ε_r = 4.4`) fills the lower box; air fills the rest.

use num_complex::Complex64;
use std::collections::HashMap;
use yee_mesh::TriMesh2D;
use yee_mom::__internal::quasi_tem_beta_for_test;

const C0: f64 = 299_792_458.0;

// ── Canonical 50Ω FR-4 microstrip (ε_r=4.4, w/h≈1.9 → Z₀≈50Ω) ──
const EPS_R: f64 = 4.4;
const H_SUB_M: f64 = 1.6e-3; // FR-4 substrate height
const W_OVER_H: f64 = 1.9; // strip width / substrate height (≈50Ω)
const FREQ_HZ: f64 = 2.0e9; // low f: box waveguide modes stay below cutoff,
// so the (no-cutoff) quasi-TEM is the dominant propagating mode.

/// Shield box width as a multiple of the strip width (a few strip widths so
/// the side walls are well clear of the strip).
const BOX_W_OVER_STRIP: f64 = 8.0;
/// Shield box height as a multiple of the substrate height (several h so the
/// top wall does not crowd the quasi-TEM field above the strip).
const BOX_H_OVER_SUB: f64 = 7.0;
/// Loose tolerance: the shielding box perturbs the open-line HJ value; a
/// larger box tightens it (CLAUDE.md placeholder-tolerance policy, ≤5–10%).
const EPS_EFF_TOL: f64 = 0.10;

/// Hammerstad-Jensen open-line effective permittivity (Balanis 3rd ed.
/// Eq. 14-1, `w/h ≥ 1` branch; identical to yee-design's Balanis
/// calculator). Inlined here because the lib-side helper is `pub(crate)`;
/// the lib unit test `eigensolver::reference::tests::
/// microstrip_hj_eps_eff_canonical_values` anchors the same closed form.
fn microstrip_eps_eff_hj(eps_r: f64, w: f64, h: f64) -> f64 {
    (eps_r + 1.0) * 0.5 + (eps_r - 1.0) * 0.5 * (1.0 + 12.0 * h / w).powf(-0.5)
}

/// Build a **shielded-microstrip** `TriMesh2D`: an `nx × ny` structured grid
/// over `[0, wb] × [0, hb]` with a rectangular **strip-conductor hole**
/// centred in x at the top of the substrate (`y ∈ [h, h+t]`). Hole cells are
/// omitted, so their border edges are mesh-boundary = PEC (the signal strip,
/// a second conductor inside the outer PEC box). Substrate (`y < h`) is
/// tagged material 1; air is tag 0.
fn shielded_microstrip_mesh(
    wb: f64,
    hb: f64,
    h_sub: f64,
    w_strip: f64,
    t_strip: f64,
    nx: usize,
    ny: usize,
) -> TriMesh2D {
    let xs: Vec<f64> = (0..=nx).map(|i| wb * (i as f64) / (nx as f64)).collect();
    let ys: Vec<f64> = (0..=ny).map(|j| hb * (j as f64) / (ny as f64)).collect();
    let xc = wb / 2.0;
    let (sx0, sx1) = (xc - w_strip / 2.0, xc + w_strip / 2.0);
    let (sy0, sy1) = (h_sub, h_sub + t_strip);
    let in_strip = |cx: f64, cy: f64| {
        cx > sx0 - 1e-12 && cx < sx1 + 1e-12 && cy > sy0 - 1e-12 && cy < sy1 + 1e-12
    };

    let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
    for &y in &ys {
        for &x in &xs {
            vertices.push([x, y]);
        }
    }
    let idx = |i: usize, j: usize| j * (nx + 1) + i;
    let mut triangles = Vec::new();
    let mut tags = Vec::new();
    for j in 0..ny {
        let yc = 0.5 * (ys[j] + ys[j + 1]);
        for i in 0..nx {
            let xcell = 0.5 * (xs[i] + xs[i + 1]);
            if in_strip(xcell, yc) {
                continue; // hole = signal-strip PEC conductor
            }
            let v00 = idx(i, j);
            let v10 = idx(i + 1, j);
            let v11 = idx(i + 1, j + 1);
            let v01 = idx(i, j + 1);
            let tag = if yc < h_sub { 1u32 } else { 0u32 };
            triangles.push([v00, v10, v11]);
            tags.push(tag);
            triangles.push([v00, v11, v01]);
            tags.push(tag);
        }
    }
    TriMesh2D::new(vertices, triangles, None, Some(tags))
        .expect("shielded microstrip mesh invariants")
}

fn fr4_eps_mu() -> (HashMap<u32, Complex64>, HashMap<u32, Complex64>) {
    let mut eps = HashMap::new();
    eps.insert(0u32, Complex64::new(1.0, 0.0)); // air
    eps.insert(1u32, Complex64::new(EPS_R, 0.0)); // FR-4
    let mut mu = HashMap::new();
    mu.insert(0u32, Complex64::new(1.0, 0.0));
    mu.insert(1u32, Complex64::new(1.0, 0.0));
    (eps, mu)
}

#[test]
fn microstrip_quasi_tem_eps_eff_matches_hammerstad_jensen() {
    // DoD-2: the quasi-TEM ε_eff = (β/k₀)² of a canonical shielded 50Ω FR-4
    // microstrip matches the Hammerstad-Jensen open-line ε_eff within a loose
    // (box-truncation-perturbed) tolerance.
    let w_strip = W_OVER_H * H_SUB_M;
    let wb = BOX_W_OVER_STRIP * w_strip;
    let hb = BOX_H_OVER_SUB * H_SUB_M;
    let (nx, ny) = (20usize, 10usize);
    let t_strip = hb / (ny as f64); // ~1 cell tall (thin signal strip)
    let mesh = shielded_microstrip_mesh(wb, hb, H_SUB_M, w_strip, t_strip, nx, ny);
    let (eps, mu) = fr4_eps_mu();

    let (beta, t_frac) = quasi_tem_beta_for_test(&mesh, &eps, &mu, FREQ_HZ)
        .expect("quasi-TEM selection must surface a transverse-dominated propagating mode");

    let k0 = std::f64::consts::TAU * FREQ_HZ / C0;
    let eps_eff_num = (beta.re / k0).powi(2);
    let eps_eff_hj = microstrip_eps_eff_hj(EPS_R, w_strip, H_SUB_M);
    let rel = (eps_eff_num - eps_eff_hj).abs() / eps_eff_hj;

    eprintln!(
        "microstrip quasi-TEM (ε_r={EPS_R}, w/h={W_OVER_H}, box {BOX_W_OVER_STRIP}w × \
         {BOX_H_OVER_SUB}h, {nx}×{ny}): β={:.4} rad/m, t-frac={t_frac:.4}, \
         ε_eff_num={eps_eff_num:.4}, ε_eff_HJ={eps_eff_hj:.4}, rel={rel:.4} (tol {EPS_EFF_TOL})",
        beta.re
    );

    // The surfaced mode must be a genuine (transverse-dominated, propagating)
    // quasi-TEM mode, not a spurious E_z / gradient capture.
    assert!(
        beta.re > 0.0 && beta.re.is_finite(),
        "quasi-TEM β must be positive-finite (got {})",
        beta.re
    );
    assert!(
        t_frac >= 0.5,
        "quasi-TEM mode must be transverse-energy-dominated (t-frac {t_frac:.4})"
    );
    // HJ reconciliation (the published-benchmark closure, loose tolerance).
    assert!(
        rel <= EPS_EFF_TOL,
        "quasi-TEM ε_eff {eps_eff_num:.4} must match Hammerstad-Jensen {eps_eff_hj:.4} within \
         {EPS_EFF_TOL} (rel {rel:.4}); a larger shield box tightens the box-truncation shift"
    );
}
