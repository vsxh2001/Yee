//! FEM-EM port-fidelity de-risk (ADR-0162 bricks **B1 + B1.5**) — a
//! power-balance + Poynting-flux energy audit on the matched straight
//! microstrip THRU.
//!
//! ## The decisive, cheap experiment
//!
//! ADR-0162's diagnosis of the filter-S21 floor (in-band peak ≈ −27 dB,
//! missing the Chebyshev mask) has a specific, testable root cause: the
//! FEM driven-sweep S-parameter extraction
//! ([`yee_fem::OpenBoundarySolver::sweep_matrix`]) normalizes the modal
//! projection by the **E-field L² self-overlap** `M = ∫|e_mode|²`, NOT by
//! the **modal power flux** `κ_m = Re ∫(e_m × h_m*)·ẑ`. The smoking gun:
//! even on the matched, lossless straight thru, the L² extraction gives
//! `|S11|² + |S21|² ≈ 0.087² + 0.778² ≈ 0.61` — ~39 % of the incident
//! power is unaccounted-for, whereas a power-conserving extraction must
//! give `≈ 1` for a lossless 2-port.
//!
//! Microstrip is **inhomogeneous** (air + dielectric), so the modal wave
//! impedance varies across the cross-section and the E-only norm
//! mis-weights the power; the dielectric region (where `√ε_r > 1`) should
//! carry more of the power flux than the air. This probe measures
//! `|S11|² + |S21|²` on the thru **two ways off the same solved FEM
//! field** — the production E-field-L² norm and a power-wave norm
//! (`κ_m = Re ∫(e_m × h_m*)·ẑ` with the quasi-TEM modal H
//! `h_t = (√ε_r/η₀)(ẑ × e_t)`) — to distinguish:
//!
//! * **normalization artifact** — power-norm lifts the balance `0.61 → ~1`
//!   ⇒ the L² norm is the magnitude bug ⇒ ADR-0162 **B2 GO**
//!   (productionize the power-norm), vs
//! * **real numerical loss** — power-norm leaves it ≈0.61 ⇒ the deficit is
//!   the known ~30 % numerical Q-floor (ADR-0156/K3) ⇒ **B2 NO-GO**, the
//!   normalization is not the magnitude fix (re-scope to the Q-floor).
//!
//! Either outcome is a valid, decisive result. This probe does **not**
//! change the production `extract_s_qp` / `extract_s11` — it adds the
//! power-norm as a separate diagnostic
//! ([`yee_fem::OpenBoundarySolver::sweep_matrix_power_balance`]) read off
//! the identical FEM field (same assembly / LU / back-substitution), so
//! the normalization is the only changed variable. B2 productionizes the
//! power-norm only if B1 confirms it here.
//!
//! ## B1.5 — the Poynting-flux audit (resolves the B1 caveat decisively)
//!
//! B1's power-norm used a quasi-TEM modal-H **approximation** and lifted
//! the balance only `0.61 → 0.67`, leaving the question ambiguous: is the
//! residual deficit real loss, or just the approximate modal H? B1.5
//! removes the approximation entirely. From the SAME solved phasor field
//! it reconstructs the **true** magnetic field `H = ∇×E / (−jωμ)` via the
//! exact Whitney-1 (Nédélec) curl `∇×N_α = 2 ∇λ_i × ∇λ_j`
//! ([`yee_fem::OpenBoundarySolver::poynting_flux_audit`]) and integrates
//! the complex Poynting vector `S = ½(E × H*)` through both port planes
//! (same 3-pt Gauss faces, outward `n̂`). The decisive ratio is
//! `P_out/P_in`: `≈1` ⇒ the solved field conserves energy ⇒ the ≈0.61
//! S-balance is an **extraction artifact** (track salvageable via a
//! flux-calibrated extraction); `≪1` ⇒ **real solve/ABC loss** (B2 NO-GO,
//! the K3 Q-floor). MEASURED: `P_out/P_in = 0.9982` — the field is
//! lossless, so the deficit is an extraction artifact (see the gate
//! docstring for the full numbers + decision).
//!
//! ## Method references (research-first)
//!
//! Modal power normalization / wave-port S-extraction: Jin, *The Finite
//! Element Method in Electromagnetics* (wave-port chapter); COMSOL RF
//! "S-Parameter Calculations" (power-flow normalization, conjugate-mode
//! overlap); arXiv 2407.21766 (`κ_m = ∫(e_m×h_m*)·ẑ`,
//! `α_i = ∫(E_tot×h_i*)·ẑ / κ_i`); Palace boundaries
//! (`S_ij = ∫E·E_inc/∫E_inc·E_inc − δ`, valid only because its wave-port
//! mode is unit-incident-**power** normalized first — the step Yee
//! skipped).
//!
//! ## GATING — CRITICAL
//!
//! Like the N2 gate this is a driven SOLVE (a per-ω sparse-LU factorization
//! on a ≲ 14 k-tet mesh plus a sub-second 2-D eigensolve per port face).
//! It is `#[ignore]`'d so the debug `cargo test --workspace` never runs it,
//! and is run only in `--release`, boxed:
//!
//! ```text
//! YEE_BOX_DIR=$(pwd) YEE_BOX_MEM=14g YEE_BOX_CPUS=3 scripts/yee-box.sh \
//!   cargo test -p yee-fem --release --test port_power_balance_001 \
//!   -- --ignored --nocapture
//! ```
//!
//! The setup mirrors the N2 straight-line gate
//! (`microstrip_eeff::fem_line_eeff_numerical_001`) EXACTLY — same
//! 6 mm × 6 mm FR-4 box on 1 mm substrate, same NX = NZ = 12 cross-section,
//! same interior-PEC trace+ground edges, same `with_coupled_whitney(true)`,
//! same production `microstrip_port_numerical` numerical-eigenmode port —
//! and uses a single short thru length (L1 = 20 mm) because the power
//! balance `|S11|² + |S21|²` is length-independent for a lossless line and
//! one solve is the cheap, decisive measurement.

#![allow(non_snake_case)]

use std::f64::consts::PI;

use nalgebra::Vector3;
use yee_fem::{
    FaceKind, MaterialDatabase, MicrostripPortGeom, OpenBoundarySolver, layered_microstrip_mesh,
    microstrip_port_numerical,
};
use yee_mesh::TetMesh3D;

// --- Geometry / material constants — identical to the N2 gate. ---
const BOX_W: f64 = 6.0e-3;
const BOX_H: f64 = 6.0e-3;
const SUB_H: f64 = 1.0e-3;
const TRACE_W: f64 = 1.0e-3;
const EPS_R: f64 = 4.4;
const NX: usize = 12;
const NZ: usize = 12;
const F_TEST: f64 = 2.0e9;

// A single short thru: L1 = 20 mm, ny = 8 (dy = 2.5 mm), the same first
// length the N2 two-length extraction uses. The power balance is
// length-independent for a lossless line, so one solve is decisive.
const L_THRU: f64 = 20.0e-3;
const NY_THRU: usize = 8;

/// Count exterior faces — same helper as the N2 gate fixture.
fn exterior_face_count(mesh: &TetMesh3D) -> usize {
    let mut face_map: std::collections::HashMap<[usize; 3], usize> =
        std::collections::HashMap::new();
    const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
    for tet in &mesh.tetrahedra {
        for &[a, b, c] in TET_FACES.iter() {
            let mut key = [tet[a], tet[b], tet[c]];
            key.sort_unstable();
            *face_map.entry(key).or_insert(0) += 1;
        }
    }
    face_map.values().filter(|&&c| c == 1).count()
}

/// Classify the two y-end-cap faces as wave-ports, everything else PEC —
/// identical to the N2 gate's `classify_microstrip_faces`.
fn classify_microstrip_faces(centroids: &[Vector3<f64>], line_len: f64) -> Vec<FaceKind> {
    let tol = 1e-9;
    centroids
        .iter()
        .map(|c| {
            if c.y < tol {
                FaceKind::WavePort(0)
            } else if (c.y - line_len).abs() < tol {
                FaceKind::WavePort(1)
            } else {
                FaceKind::Pec
            }
        })
        .collect()
}

fn numerical_port_geom() -> MicrostripPortGeom {
    MicrostripPortGeom {
        trace_w: TRACE_W,
        sub_h: SUB_H,
        eps_r: EPS_R,
        box_w: BOX_W,
        box_h: BOX_H,
    }
}

/// Build the N2 numerical-port thru solver for a line of length `line_len`
/// (ny longitudinal cells), borrowing the caller-owned `mesh` so the test
/// can run **both** the power-balance sweep (B1) and the Poynting-flux
/// audit (B1.5) on the SAME solver / SAME solved field.
///
/// Bit-identical mesh / interior-PEC / coupled-Whitney / numerical-port
/// path to `microstrip_eeff::solve_line_numerical`. The mesh + material
/// database are built by [`build_thru_mesh`] and owned by the test; this
/// returns the solver that borrows them.
fn build_thru_solver<'m>(
    mesh: &'m TetMesh3D,
    material_db: MaterialDatabase,
    line_len: f64,
    ground_pred: &dyn Fn(Vector3<f64>, Vector3<f64>) -> bool,
    trace_pred: &dyn Fn(Vector3<f64>, Vector3<f64>) -> bool,
    geom: &MicrostripPortGeom,
) -> OpenBoundarySolver<'m> {
    let n_exterior = exterior_face_count(mesh);
    let picker = OpenBoundarySolver::new(
        mesh,
        vec![FaceKind::Pec; n_exterior],
        Vec::new(),
        MaterialDatabase::new(),
    )
    .expect("picker solver must build");
    let ground_edges = picker.interior_edges_matching(ground_pred);
    let trace_edges = picker.interior_edges_matching(trace_pred);
    let mut interior_pec: Vec<usize> = ground_edges;
    interior_pec.extend(trace_edges.iter().copied());
    interior_pec.sort_unstable();
    interior_pec.dedup();
    assert!(
        !trace_edges.is_empty(),
        "trace_pred must select at least one interior edge on the z = sub_h trace footprint"
    );
    let centroids = picker.exterior_face_centroids();
    let kinds = classify_microstrip_faces(&centroids, line_len);
    drop(picker);

    let port_0 = microstrip_port_numerical(geom, F_TEST)
        .expect("numerical-eigenmode port (face 0) must build");
    let port_1 = microstrip_port_numerical(geom, F_TEST)
        .expect("numerical-eigenmode port (face 1) must build");
    OpenBoundarySolver::new(mesh, kinds, vec![port_0, port_1], material_db)
        .expect("two-port numerical-port microstrip solver must build")
        .with_interior_pec_edges(interior_pec.iter().copied())
        .with_coupled_whitney(true)
}

/// Per-point relative permittivity for the quasi-TEM modal-H relation
/// `h_t = (√ε_r/η₀)(ẑ × e_t)`: the layered microstrip is FR-4 in the
/// substrate (z < sub_h) and air above. The dielectric region weights
/// √ε_r heavier — the spatial weighting the homogeneous E-only L² norm
/// misses. On the y-end-cap port face the transverse plane is (x, z), so
/// z is the substrate-normal axis.
fn eps_r_at(p: Vector3<f64>) -> f64 {
    if p.z < SUB_H { EPS_R } else { 1.0 }
}

/// **THE B1 DE-RISK PROBE** (ADR-0162). Drives the matched straight FR-4
/// microstrip thru through the production numerical-eigenmode wave-port
/// and prints the matched-thru power balance `|S11|² + |S21|²` computed
/// **two ways off the identical FEM field**:
///
/// 1. the production **E-field-L² normalization** (`s_l2`, the
///    smoking-gun ≈0.61 ADR-0162 reports), and
/// 2. a **power-wave normalization** (`s_power`,
///    `κ_m = Re∫(e_m×h_m*)·ẑ` with the quasi-TEM modal H).
///
/// **B1.5 (the decisive caveat-resolver).** The B1 power-norm used a
/// quasi-TEM modal-H *approximation* and lifted the balance only
/// `0.61 → 0.67`, leaving a caveat: is the residual deficit real loss, or
/// an artifact of that approximate modal H? B1.5 removes the approximation
/// — it reconstructs the **true** magnetic field `H = ∇×E/(−jωμ)` from the
/// solved electric DoFs (exact Whitney-1 curl, `poynting_flux_audit`) and
/// integrates the complex Poynting vector `S = ½(E×H*)` through the port
/// faces. The decisive ratio `P_out/P_in`:
///
/// * `≈ 1` ⇒ the **solved field conserves energy** ⇒ the ≈0.61 S-balance
///   is an **extraction artifact** ⇒ port-fidelity track SALVAGEABLE
///   (flux-calibrated extraction recovers it); vs
/// * `≪ 1` ⇒ **real volume/ABC loss** in the solve ⇒ B2 NO-GO confirmed,
///   the floor is numerical/ABC dissipation (the K3 Q-floor).
///
/// `P_out/P_in` is robust to the H sign/`(ωμ)`-scaling convention (both
/// fluxes use the same H, so the convention cancels in the ratio) and to
/// any global field normalization — only a *relative* error between the
/// two identically-reconstructed port faces could bias it.
///
/// The probe does NOT assert a target balance/ratio — it is a measurement
/// that drives the decision; it only asserts the pipeline did not
/// degenerate (finite outputs, non-zero |S21| / P_in), so a hard-wall /
/// collapsed-port regression surfaces rather than silently reporting a
/// meaningless number.
///
/// ## MEASURED RESULT (boxed --release, base a64fe86)
///
/// ```text
///   B1   normalization      |S11|    |S21|    |S11|²+|S21|²
///        E-field L²         0.0956   0.7780   0.6145   <- ADR-0162 smoking gun (≈0.61) CONFIRMED
///        power-wave (√ε_r)  0.0646   0.8175   0.6725   <- quasi-TEM-h approx: barely lifts (+0.058)
///
///   B1.5 Poynting-flux energy audit (TRUE H = ∇×E/(−jωμ), no modal approx)
///        P_in  (into driven port 0)  = 1.7137e-10 W
///        P_out (out of port 1)       = 1.7107e-10 W
///        P_out / P_in  (FIELD ratio) = 0.9982   <- THE DECISIVE NUMBER
///        |S21|²/(1−|S11|²) (S ratio) = 0.6109   <- the SAME field, via the S-params
/// ```
///
/// **Decision: EXTRACTION ARTIFACT — the port-fidelity track is
/// SALVAGEABLE.** The B1.5 audit is decisive and overturns the B1-only
/// reading. The solved field's true Poynting power balance is
/// `P_out/P_in = 0.9982` — it transmits essentially **all** incident power
/// from port 0 to port 1, losing ~0.2 %. So the ≈0.61 S-parameter balance
/// is **NOT real loss**: the field conserves energy; the S-parameter
/// *extraction* under-counts the transmitted power. The proof is the
/// side-by-side `|S21|²/(1−|S11|²) = 0.6109` (the same field through the
/// S-formula) vs `P_out/P_in = 0.9982` (the same field through Poynting) —
/// a ~0.39 gap that is purely the extraction normalization.
///
/// This resolves the B1 caveat (B1's quasi-TEM power-norm lifted the
/// balance only +0.06, leaving it ambiguous whether the residual was real
/// loss or an approximate-modal-H artifact): with the **true** H there is
/// **no** loss, so the residual is an extraction artifact. The B1 √ε_r
/// power-norm helped little only because its quasi-TEM modal-H *shape* is
/// itself approximate, not because the deficit is loss. The fix is a
/// **flux-calibrated extraction** (normalize the modal projection so it
/// counts the true modal power flux `½Re∫(e_m×h_m*)·ẑ` consistently) —
/// reconsider B2 in that form rather than as the quasi-TEM √ε_r reweight.
///
/// `P_out/P_in` is robust to the H sign / `(ωμ)`-scaling convention (it
/// cancels in the ratio) and to global field normalization, so the
/// conclusion does not hinge on a convention choice.
#[test]
#[ignore = "driven SOLVE (per-ω sparse LU + 2-D eigensolves); run only in --release, boxed"]
fn port_power_balance_001() {
    let geom = numerical_port_geom();
    let omega = 2.0 * PI * F_TEST;

    // Build the mesh ONCE and keep it alive so the solver (which borrows
    // it) can run BOTH the B1 power-balance sweep and the B1.5
    // Poynting-flux audit on the SAME solver / SAME solved field.
    let (mesh, material_db, ground_pred, trace_pred) =
        layered_microstrip_mesh(BOX_W, BOX_H, L_THRU, SUB_H, TRACE_W, NX, NY_THRU, NZ)
            .expect("layered_microstrip_mesh must build for the chosen geometry");
    let solver = build_thru_solver(&mesh, material_db, L_THRU, &ground_pred, &trace_pred, &geom);

    // ── B1: two S-matrices (L² + power-wave) off one solve. ──
    let sweep = solver
        .sweep_matrix_power_balance(&[omega], &eps_r_at)
        .expect("power-balance diagnostic sweep must succeed");

    // ── B1.5: Poynting-flux energy audit of the driven (port 0) solve,
    // using the TRUE H = ∇×E/(−jωμ) — no modal-H approximation. ──
    let audit = solver
        .poynting_flux_audit(omega, 0)
        .expect("Poynting-flux audit must succeed");

    // Single frequency, 2-port: read S11 = s[(0,0)], S21 = s[(1,0)].
    let s_l2 = &sweep.s_l2[0];
    let s_pow = &sweep.s_power[0];

    let s11_l2 = s_l2[(0, 0)];
    let s21_l2 = s_l2[(1, 0)];
    let s11_pow = s_pow[(0, 0)];
    let s21_pow = s_pow[(1, 0)];

    let bal_l2 = s11_l2.norm_sqr() + s21_l2.norm_sqr();
    let bal_pow = s11_pow.norm_sqr() + s21_pow.norm_sqr();

    // The field's own power-transmission ratio, for cross-check: the
    // S-parameter analogue of the Poynting ratio is |S21|²/(1−|S11|²).
    let s_param_ratio = s21_l2.norm_sqr() / (1.0 - s11_l2.norm_sqr());

    eprintln!(
        "\n==== B1 + B1.5 PORT-FIDELITY PROBE (ADR-0162) ====\n\
         thru: L = {:.1} mm, 6×6 mm FR-4 box, numerical-eigenmode port, f = {:.2} GHz\n\
         \n\
         --- B1: E-field L² normalization (production extract_s_qp) ---\n\
         |S11|_L2          : {:.4}\n\
         |S21|_L2          : {:.4}\n\
         |S11|²+|S21|²_L2  : {bal_l2:.4}   <-- the ADR-0162 smoking gun (≈0.61?)\n\
         \n\
         --- B1: power-wave normalization (κ_m = Re∫(e_m×h_m*)·ẑ, quasi-TEM h) ---\n\
         |S11|_pow         : {:.4}\n\
         |S21|_pow         : {:.4}\n\
         |S11|²+|S21|²_pow : {bal_pow:.4}   (quasi-TEM-h approx; lifts only ~+0.06)\n\
         \n\
         --- B1.5: Poynting-flux energy audit (TRUE H = ∇×E/(−jωμ), no modal approx) ---\n\
         P_in  (into driven port 0)  : {:.6e} W\n\
         P_out (out of port 1)       : {:.6e} W\n\
         per-port net leaving (W)    : [{:.4e}, {:.4e}]\n\
         P_out/P_in  (FIELD ratio)   : {:.4}   <-- THE DECISIVE NUMBER (→1 or ≪1?)\n\
         |S21|²/(1−|S11|²) (S ratio) : {s_param_ratio:.4}   (S-parameter analogue, cross-check)\n\
         \n\
         decision: {}\n\
         ==================================================",
        L_THRU * 1e3,
        F_TEST / 1e9,
        s11_l2.norm(),
        s21_l2.norm(),
        s11_pow.norm(),
        s21_pow.norm(),
        audit.p_in,
        audit.p_out,
        audit.p_leaving[0],
        audit.p_leaving[1],
        audit.power_ratio,
        if audit.power_ratio > 0.95 {
            "P_out/P_in ≈ 1 ⇒ the SOLVED FIELD conserves energy ⇒ the ≈0.61 S-balance is an \
             EXTRACTION ARTIFACT ⇒ port-fidelity track SALVAGEABLE (flux-calibrated extraction)"
        } else if audit.power_ratio < 0.8 {
            "P_out/P_in ≪ 1 ⇒ REAL volume/ABC loss in the solve ⇒ the floor is numerical/ABC \
             dissipation (K3 Q-floor) ⇒ B2 NO-GO confirmed, re-scope to the solver/ABC"
        } else {
            "P_out/P_in in (0.8, 0.95) ⇒ PARTIAL solve loss ⇒ mixed cause; inspect before B2"
        },
    );

    // Non-degeneracy only — this is a MEASUREMENT probe. It does NOT assert
    // a target balance or ratio (that would pre-judge the GO/NO-GO the
    // numbers must drive). Surface a pipeline collapse, then let the
    // printed numbers speak.
    assert!(
        bal_l2.is_finite() && bal_pow.is_finite() && audit.power_ratio.is_finite(),
        "a probe output is non-finite (bal_L²={bal_l2}, bal_pow={bal_pow}, \
         P_out/P_in={}) — the diagnostic diverged and cannot conclude",
        audit.power_ratio
    );
    assert!(
        s21_l2.norm() > 1e-6 && s21_pow.norm() > 1e-6 && audit.p_in.abs() > 0.0,
        "|S21| or P_in collapsed to ~0 (|S21|_L²={:.3e}, P_in={:.3e} W) — the numerical \
         port degenerated to a hard wall; re-run/inspect before reading the ratios",
        s21_l2.norm(),
        audit.p_in,
    );
}
