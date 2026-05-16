//! Port excitation, dense LU, and S-parameter extraction.

#![allow(dead_code)]
// PlanarMoM::run (Task 10) consumes these.

use crate::basis::RwgBasis;
use crate::fill::impedance_matrix;
use crate::greens::FreeSpaceGreen;
use crate::iterative::{GmresParams, gmres_jacobi};
use crate::ports::{DeltaGapPort, Port};
use crate::roughness::{RoughnessModel, SIGMA_COPPER};
use faer::Mat;
use faer::linalg::solvers::{PartialPivLu, Solve};
use num_complex::Complex64;
use yee_core::{Error, FreqRange};
use yee_io::touchstone::{File as TouchstoneFile, Format, FreqUnit};

/// Selector for the linear solver used by [`s_parameters_at_freq_with_solver`].
///
/// The default ([`LinearSolver::Direct`]) is the dense partial-pivot LU
/// path used since Phase 1.0 — appropriate for `n ≲ 10⁴`. For large
/// problems (`n ≥ 5 × 10⁴` is the rough crossover where dense LU exceeds
/// 32 GB of GPU memory) use [`LinearSolver::GmresJacobi`].
#[derive(Debug, Clone, Copy, Default)]
pub enum LinearSolver {
    /// Dense partial-pivot LU via `faer::linalg::solvers::PartialPivLu`.
    #[default]
    Direct,
    /// Restarted GMRES(m) with diagonal Jacobi preconditioning. See
    /// [`crate::iterative::gmres_jacobi`].
    GmresJacobi(GmresParams),
}

/// Apply a frequency-dependent surface-roughness loss multiplier to the
/// impedance matrix in-place.
///
/// Phase 1.4.0 walking skeleton: scales every entry of `z` by the same
/// `K(f)` returned by `roughness.loss_multiplier(freq_hz, SIGMA_COPPER)`,
/// which is mathematically equivalent to scaling the final `Z_in` (and
/// therefore `S11` through the bilinear transform) by `K(f)`. This is a
/// coarse approximation — see the limitation documented in
/// [`crate::roughness`]. Conductivity is hard-coded to copper here;
/// Phase 1.4.1 will plumb a material-aware sigma through the basis.
pub(crate) fn apply_roughness(z: &mut Mat<Complex64>, roughness: &RoughnessModel, freq_hz: f64) {
    let k = roughness.loss_multiplier(freq_hz, SIGMA_COPPER);
    let scale = Complex64::new(k, 0.0);
    let (nrows, ncols) = (z.nrows(), z.ncols());
    for j in 0..ncols {
        for i in 0..nrows {
            z[(i, j)] *= scale;
        }
    }
}

/// Build a 1 V delta-gap RHS across every edge tagged with `port_tag`.
/// `b[k] = V × length_k` for port edges, zero elsewhere.
///
/// Thin compatibility wrapper around [`DeltaGapPort`] with `voltage = 1 + 0i`.
/// Retained so existing internal/test call sites that pre-date the Phase 1.3
/// [`Port`] trait keep compiling against a stable signature.
pub(crate) fn delta_gap_rhs(basis: &RwgBasis, port_tag: u32) -> Mat<Complex64> {
    let port = DeltaGapPort {
        tag: port_tag,
        voltage: Complex64::new(1.0, 0.0),
    };
    // Free-space delta-gap RHS is frequency-independent; the `freq_hz`
    // argument exists only for ports whose RHS depends on it (wave ports
    // in Phase 1.3.1+). Passing 0.0 here is safe.
    port.rhs(basis, 0.0)
}

/// Solve at a single frequency and return S11 referenced to `z0_ref`.
///
/// Pipeline:
/// 1. Build the free-space Green's function at `freq_hz`.
/// 2. Assemble the MPIE impedance matrix `Z` over the RWG basis.
/// 3. Build the RHS `b = port.rhs(basis, freq_hz)`.
/// 4. Factorise `Z = L U` (partial pivoting) and solve `Z i = b`.
/// 5. Recover the port current as `port.port_current(basis, i)`. For a
///    delta-gap or uniform-mode wave port this is the Galerkin projection
///    `Σ_k length_k · i_k` over port edges, matching the inner product
///    that defines `b`.
/// 6. Return `S11 = (Z_in − Z0) / (Z_in + Z0)` with `Z_in = V / I_port` and
///    `V = port.port_voltage()`.
///
/// Returns [`Error::Numerical`] if `|I_port|` drops below `1e-30`, which
/// indicates a pathological port tagging (no real current driven).
pub(crate) fn s_parameters_at_freq(
    basis: &RwgBasis,
    port: &dyn Port,
    freq_hz: f64,
    z0_ref: f64,
    roughness: Option<&RoughnessModel>,
) -> Result<Complex64, Error> {
    s_parameters_at_freq_with_solver(
        basis,
        port,
        freq_hz,
        z0_ref,
        roughness,
        LinearSolver::Direct,
    )
}

/// Same numerics as [`s_parameters_at_freq`] but with an explicit solver
/// selection. The direct LU path is the canonical Phase 1.0 numerics and
/// remains the default for callers that don't care; the GMRES path is
/// intended for `n ≥ 50k` where dense LU overflows GPU memory.
///
/// Returns [`Error::Numerical`] for both the pathological-port case and
/// for GMRES non-convergence (the iterative path surfaces the final
/// relative residual in the error message).
pub(crate) fn s_parameters_at_freq_with_solver(
    basis: &RwgBasis,
    port: &dyn Port,
    freq_hz: f64,
    z0_ref: f64,
    roughness: Option<&RoughnessModel>,
    solver: LinearSolver,
) -> Result<Complex64, Error> {
    let green = FreeSpaceGreen::new(freq_hz);
    let mut z = impedance_matrix(basis, &green);
    if let Some(rough) = roughness {
        apply_roughness(&mut z, rough, freq_hz);
    }
    let b = port.rhs(basis, freq_hz);

    let i = match solver {
        LinearSolver::Direct => {
            // Partial-pivot LU is the canonical dense solver for the
            // (non-Hermitian but symmetric) MPIE impedance matrix; full
            // pivoting is overkill and QR would needlessly inflate the
            // constant factor on the O(N^3) work.
            let lu = PartialPivLu::new(z.as_ref());
            lu.solve(b.as_ref())
        }
        LinearSolver::GmresJacobi(params) => {
            let res = gmres_jacobi(z.as_ref(), b.as_ref(), params);
            if !res.converged {
                return Err(Error::Numerical(format!(
                    "GMRES failed to converge: final residual {:.3e} after {} iterations",
                    res.final_residual, res.iterations
                )));
            }
            res.x
        }
    };

    let i_port = port.port_current(basis, &i);

    if i_port.norm() < 1e-30 {
        return Err(Error::Numerical(
            "port current vanished; check port tagging".into(),
        ));
    }

    let v_port = port.port_voltage();
    let z_in = v_port / i_port;
    let z0 = Complex64::new(z0_ref, 0.0);
    Ok((z_in - z0) / (z_in + z0))
}

/// Run the sweep and produce a Touchstone file (Format::RealImag, FreqUnit::Hz).
///
/// One `s_parameters_at_freq` call per frequency in `freq_range`. The
/// resulting `TouchstoneFile` is filled with all fields populated — z0,
/// freq_unit, format, n_ports, freq_hz, data, and a single attribution
/// comment — so it round-trips through `yee_io::touchstone::write` without
/// further massaging.
pub(crate) fn s_parameters_sweep(
    basis: &RwgBasis,
    port: &dyn Port,
    freq_range: FreqRange,
    z0_ref: f64,
    roughness: Option<&RoughnessModel>,
) -> Result<TouchstoneFile, Error> {
    let mut freq_hz = Vec::new();
    let mut data: Vec<Vec<Complex64>> = Vec::new();
    for f in freq_range.iter() {
        let s11 = s_parameters_at_freq(basis, port, f, z0_ref, roughness)?;
        freq_hz.push(f);
        data.push(vec![s11]);
    }
    Ok(TouchstoneFile {
        n_ports: 1,
        z0: z0_ref,
        freq_unit: FreqUnit::Hz,
        format: Format::RealImag,
        freq_hz,
        data,
        comments: vec!["! Generated by yee-mom Phase 1.0 free-space dipole solver".to_string()],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::RwgBasis;
    use nalgebra::Vector3;
    use yee_mesh::TriMesh;

    fn two_tri_mesh_with_port() -> TriMesh {
        let vertices = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.1, 0.0, 0.0),
            Vector3::new(0.1, 0.1, 0.0),
            Vector3::new(0.0, 0.1, 0.0),
        ];
        let triangles = vec![[0u32, 1, 2], [0u32, 2, 3]];
        // Different non-zero tags → shared diagonal edge is the port edge
        // under the "boundary between two tagged regions" port rule.
        // `port_basis_indices(1)` matches because `min(1, 2) == 1`.
        let tags = vec![1u32, 2u32];
        TriMesh::new(vertices, triangles, tags).unwrap()
    }

    #[test]
    fn delta_gap_rhs_length_weighting() {
        let basis = RwgBasis::from_mesh(two_tri_mesh_with_port()).unwrap();
        let b = delta_gap_rhs(&basis, 1);
        let port_indices: Vec<usize> = basis.port_basis_indices(1).collect();
        assert!(!port_indices.is_empty(), "expected at least one port edge");
        for k in port_indices {
            let expected = basis.edges[k].length;
            assert!((b[(k, 0)].re - expected).abs() < 1e-12);
            assert!(b[(k, 0)].im.abs() < 1e-12);
        }
    }

    #[test]
    fn sweep_produces_n_points_rows() {
        let basis = RwgBasis::from_mesh(two_tri_mesh_with_port()).unwrap();
        let freq = FreqRange::new(1.0e9, 1.5e9, 3).unwrap();
        let port = DeltaGapPort {
            tag: 1,
            voltage: Complex64::new(1.0, 0.0),
        };
        let file = s_parameters_sweep(&basis, &port, freq, 50.0, None).expect("sweep");
        assert_eq!(file.freq_hz.len(), 3);
        assert_eq!(file.data.len(), 3);
        assert_eq!(file.n_ports, 1);
        assert!(matches!(file.format, Format::RealImag));
        assert!(matches!(file.freq_unit, FreqUnit::Hz));
    }
}
