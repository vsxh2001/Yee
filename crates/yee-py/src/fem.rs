//! `yee.fem` submodule — 3-D FEM eigenmode binding.
//!
//! Python wrapper for the Phase 4.fem.eig.0 walking-skeleton pipeline:
//! a one-shot `solve_cavity(a, b, d, nx, ny, nz, num_eigs)` that builds a
//! Kuhn 6-tet uniform mesh of a closed PEC rectangular cavity via
//! [`yee_mesh::TetMesh3D::cavity_uniform`], assembles the Nedelec
//! curl-curl pencil via [`yee_fem::FemEigenAssembly::new_free_space`],
//! and solves `K e = k² M e` for the lowest `num_eigs` physical modes via
//! shift-invert deflated inverse-power iteration
//! ([`yee_fem::InverseIterEigen`]).
//!
//! Mirrors the Rust-side `yee_validation::run_fem_eig_001_rectangular_cavity`
//! driver (Phase 4 T7) end-to-end; the shift heuristic is the same
//! `σ = 2.5 · k₀_TE101²` value vetted on the (8, 6, 10) mesh in T7's
//! production-gate run. See the documented limitation in the validation
//! driver's source for why this shift sits above the lowest physical
//! mode; the [`yee_fem::SparseEigen`] trait keeps the choice behind an
//! abstraction so a future LOBPCG / ARPACK swap fixes the dependency in
//! one PR.
//!
//! ## Submodule registration
//!
//! Mirrors the `yee.touchstone` / `yee.eigensolver` pattern: we insert
//! the module into `sys.modules` from `lib.rs` so that
//! `from yee.fem import solve_cavity` works in addition to the
//! attribute-access form `yee.fem.solve_cavity`.

use num_complex::Complex64;
use numpy::{IntoPyArray, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyComplex, PyDict, PyList};
use std::f64::consts::PI;
use yee_core::units::C0;
use yee_fem::{
    DispersiveError, DispersiveSolver, FemEigenAssembly, InverseIterEigen, Material,
    MaterialDatabase, MaterialPole, SparseEigen,
};
use yee_mesh::{MaterialTag, TetMesh3D};

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(solve_cavity, m)?)?;
    m.add_function(wrap_pyfunction!(solve_cavity_dispersive, m)?)?;
    Ok(())
}

/// Solve the lowest-`num_eigs` resonant modes of a closed PEC
/// rectangular cavity.
///
/// Parameters
/// ----------
/// a, b, d : float
///     Cavity extents along the `x`, `y`, `z` axes, in metres. Must be
///     strictly positive and finite.
/// nx, ny, nz : int
///     Brick subdivisions along the `x`, `y`, `z` axes. Each axis-aligned
///     brick is decomposed into 6 tetrahedra (Kuhn decomposition), so the
///     total tet count is `nx * ny * nz * 6`. Must be ≥ 1.
/// num_eigs : int, optional
///     Number of lowest physical eigenvalues to return. Default 10.
///
/// Returns
/// -------
/// frequencies : list[float]
///     Resonant frequencies in Hz, length `num_eigs`, sorted ascending.
/// mode_coefficients : numpy.ndarray
///     `(num_interior_edges, num_eigs)` `float64` array of eigenvectors
///     on the interior-edge basis (PEC-eliminated). Column `n` is the
///     mode for `frequencies[n]`.
///
/// Notes
/// -----
/// The shift `σ = 2.5 · k₀_TE101²` is computed from the analytic
/// `TE_{101}` wavenumber `k₀_TE101² = (π/a)² + (π/d)²`, matching the
/// Phase 4 T7 production-gate driver. The current
/// `yee_fem::InverseIterEigen` is the Phase 4 T5 escape-hatch
/// implementation; a future LOBPCG / ARPACK swap on the
/// `yee_fem::SparseEigen` trait will remove the shift-sensitivity caveat
/// documented in `yee_validation::run_fem_eig_001_rectangular_cavity`.
///
/// Raises
/// ------
/// ValueError
///     If any extent is non-positive / non-finite, any subdivision is
///     zero, or `num_eigs == 0`.
/// RuntimeError
///     If the inner sparse LU or any eigenmode fails to converge in the
///     solver's iteration budget.
#[pyfunction]
#[pyo3(signature = (a, b, d, nx, ny, nz, num_eigs = 10))]
#[allow(clippy::too_many_arguments)]
fn solve_cavity(
    py: Python<'_>,
    a: f64,
    b: f64,
    d: f64,
    nx: usize,
    ny: usize,
    nz: usize,
    num_eigs: usize,
) -> PyResult<(Vec<f64>, Py<PyArray2<f64>>)> {
    // ---- 1. Build the cavity mesh -----------------------------------
    let mesh =
        TetMesh3D::cavity_uniform(a, b, d, nx, ny, nz).map_err(crate::errors::yee_mesh_to_py)?;

    // ---- 2. Assemble free-space K, M with PEC Dirichlet -------------
    let assembled = FemEigenAssembly::new_free_space(&mesh)
        .assemble()
        .map_err(crate::errors::yee_to_py)?;

    // ---- 3. Shift heuristic — same as the Phase 4 T7 driver ---------
    //
    // `solve_cavity` cannot know which returned eigenvalue is TE_{101}
    // ahead of time, so the shift uses the analytic TE_{101} wavenumber
    // for the generic rectangular cavity: `k₀² = (π/a)² + (π/d)²` (the
    // `(m, n, p) = (1, 0, 1)` mode is insensitive to `b`). The factor
    // `2.5` lifts `σ` above the gradient-kernel cluster at `k² = 0`
    // and above the lowest physical mode, per the documented limitation
    // of the inverse-power escape-hatch solver — see
    // `yee_validation::run_fem_eig_001_rectangular_cavity` source.
    let pi = std::f64::consts::PI;
    let k0_te101_sq = (pi / a).powi(2) + (pi / d).powi(2);
    let sigma = 2.5 * k0_te101_sq;

    // ---- 4. Solve K e = k² M e --------------------------------------
    let pairs = InverseIterEigen::default()
        .solve(&assembled.k, &assembled.m, num_eigs, sigma)
        .map_err(crate::errors::yee_to_py)?;

    // ---- 5. Convert k² → frequency, sort with paired eigenvectors ---
    //
    // `pairs.k` is sorted ascending by `k²` already (see
    // `InverseIterEigen::solve` post-processing), but we sort defensively
    // here so the returned `frequencies` list contract is
    // implementation-agnostic. The eigenvectors are reordered in
    // lockstep so column `n` of `mode_coefficients` always matches
    // `frequencies[n]`.
    let n_modes = pairs.k.len();
    let n_dofs = pairs.e.nrows();
    let mut order: Vec<usize> = (0..n_modes).collect();
    order.sort_by(|&i, &j| pairs.k[i].total_cmp(&pairs.k[j]));

    let two_pi = 2.0 * pi;
    let frequencies: Vec<f64> = order
        .iter()
        .map(|&i| {
            let k_sq = pairs.k[i];
            let k_abs = if k_sq > 0.0 { k_sq.sqrt() } else { 0.0 };
            C0 * k_abs / two_pi
        })
        .collect();

    // Build a `(n_dofs, n_modes)` row-major buffer in the permuted
    // column order.
    let mut buf: Vec<f64> = Vec::with_capacity(n_dofs * n_modes);
    for row in 0..n_dofs {
        for &col in &order {
            buf.push(pairs.e[(row, col)]);
        }
    }
    let modes = numpy::ndarray::Array2::from_shape_vec((n_dofs, n_modes), buf)
        .expect("buffer length matches (n_dofs, n_modes) by construction")
        .into_pyarray(py)
        .unbind();

    Ok((frequencies, modes))
}

/// Extract a single `MaterialPole` from a Python pole dict.
///
/// Accepted shapes:
///
/// * `{"kind": "drude",   "omega_p": float, "gamma": float}`
/// * `{"kind": "lorentz", "omega_0": float, "omega_p": float, "gamma": float}`
/// * `{"kind": "debye",   "eps_s_minus_eps_inf": float, "tau": float}`
///
/// The `kind` string is case-insensitive. Missing or wrong-typed fields
/// surface as a `ValueError` with the offending key.
fn pole_from_pydict(pole_dict: &Bound<'_, PyDict>) -> PyResult<MaterialPole> {
    let kind: String = pole_dict
        .get_item("kind")?
        .ok_or_else(|| PyValueError::new_err("pole dict missing required key 'kind'"))?
        .extract()
        .map_err(|e| PyValueError::new_err(format!("pole 'kind' must be a string: {e}")))?;
    let kind_lc = kind.to_ascii_lowercase();

    let get_f64 = |key: &str| -> PyResult<f64> {
        let val = pole_dict
            .get_item(key)?
            .ok_or_else(|| PyValueError::new_err(format!("pole '{kind}' missing key '{key}'")))?;
        val.extract::<f64>().map_err(|e| {
            PyValueError::new_err(format!("pole '{kind}' key '{key}' must be float: {e}"))
        })
    };

    match kind_lc.as_str() {
        "drude" => Ok(MaterialPole::Drude {
            omega_p: get_f64("omega_p")?,
            gamma: get_f64("gamma")?,
        }),
        "lorentz" => Ok(MaterialPole::Lorentz {
            omega_0: get_f64("omega_0")?,
            omega_p: get_f64("omega_p")?,
            gamma: get_f64("gamma")?,
        }),
        "debye" => Ok(MaterialPole::Debye {
            eps_s_minus_eps_inf: get_f64("eps_s_minus_eps_inf")?,
            tau: get_f64("tau")?,
        }),
        other => Err(PyValueError::new_err(format!(
            "unknown pole kind '{other}': expected one of 'drude', 'lorentz', 'debye'"
        ))),
    }
}

/// Extract a single `(MaterialTag, Material)` pair from a Python material
/// dict shaped as `{"tag": int, "eps_inf": float, "mu_r": float, "poles": list[dict]}`.
///
/// `poles` is optional and defaults to an empty list (non-dispersive
/// constant-`ε_∞` material).
fn material_from_pydict(mat_dict: &Bound<'_, PyDict>) -> PyResult<(MaterialTag, Material)> {
    let tag: MaterialTag = mat_dict
        .get_item("tag")?
        .ok_or_else(|| PyValueError::new_err("material dict missing required key 'tag'"))?
        .extract()
        .map_err(|e| PyValueError::new_err(format!("material 'tag' must be int: {e}")))?;
    let eps_inf: f64 = mat_dict
        .get_item("eps_inf")?
        .ok_or_else(|| PyValueError::new_err("material dict missing required key 'eps_inf'"))?
        .extract()
        .map_err(|e| PyValueError::new_err(format!("material 'eps_inf' must be float: {e}")))?;
    let mu_r: f64 = mat_dict
        .get_item("mu_r")?
        .ok_or_else(|| PyValueError::new_err("material dict missing required key 'mu_r'"))?
        .extract()
        .map_err(|e| PyValueError::new_err(format!("material 'mu_r' must be float: {e}")))?;

    let poles: Vec<MaterialPole> = match mat_dict.get_item("poles")? {
        None => Vec::new(),
        Some(py_poles) => {
            let list = py_poles.cast::<PyList>().map_err(|e| {
                PyValueError::new_err(format!("material 'poles' must be a list of dicts: {e}"))
            })?;
            let mut out = Vec::with_capacity(list.len());
            for item in list.iter() {
                let pole_dict = item.cast::<PyDict>().map_err(|e| {
                    PyValueError::new_err(format!("material 'poles' entry must be a dict: {e}"))
                })?;
                out.push(pole_from_pydict(pole_dict)?);
            }
            out
        }
    };

    Ok((
        tag,
        Material {
            eps_inf,
            mu_r,
            poles,
        },
    ))
}

/// Solve a dispersive complex eigenfrequency for a closed PEC rectangular
/// cavity filled with one or more dispersive materials.
///
/// Wraps [`yee_fem::DispersiveSolver::solve_with_newton`] (Phase 4.fem.eig.1
/// D5) into a Python-ergonomic dict-in / dict-out shape. Builds a Kuhn 6-tet
/// uniform mesh via [`yee_mesh::TetMesh3D::cavity_uniform`], constructs a
/// [`yee_fem::MaterialDatabase`] from the Python `materials` list, and runs
/// the outer fixed-point Newton ω-tracker from `omega_warm_start_hz`.
///
/// Parameters
/// ----------
/// a, b, d : float
///     Cavity extents along the `x`, `y`, `z` axes, in metres. Must be
///     strictly positive and finite.
/// nx, ny, nz : int
///     Brick subdivisions along the `x`, `y`, `z` axes (Kuhn 6-tet
///     decomposition; total tet count is `nx * ny * nz * 6`). Must be ≥ 1.
/// materials : list[dict]
///     Per-tag dispersive materials. Each entry must be a dict with keys
///     ``tag`` (int — MaterialTag, matches
///     `mesh.tetrahedron_material[i]`), ``eps_inf`` (float — high-frequency
///     ε), ``mu_r`` (float — real μ_r), and an optional ``poles`` (list of
///     dicts). For `cavity_uniform` every tet currently carries the bulk
///     tag ``0`` (see `yee_fem::dispersive::BULK_TAG`), so the standard
///     input is ``materials = [{"tag": 0, "eps_inf": ..., "mu_r": ...,
///     "poles": [...]}]``. Each pole dict carries ``kind`` (one of
///     ``"drude"``, ``"lorentz"``, ``"debye"``) plus the per-kind
///     parameters — see :func:`yee.fem` module docs.
/// omega_warm_start_hz : float
///     Real-valued warm-start frequency in Hz. Internally converted to
///     `ω₀ = 2π · f + 0j` (rad/s) for the Newton tracker. Typical value:
///     the air-resonance frequency from
///     :func:`yee.fem.solve_cavity` on the same geometry.
/// max_iter : int, optional
///     Inner-solver per-mode iteration cap forwarded to
///     [`yee_fem::ComplexInverseIterEigen`]. Default ``8`` (also caps the
///     outer Newton budget; see Notes).
/// tol : float, optional
///     Inner-solver Rayleigh-quotient tolerance. Default ``1e-6``.
///
/// Returns
/// -------
/// dict
///     Mapping with the keys:
///
///     * ``"frequency_hz"`` (:class:`complex`) — converged complex
///       eigenfrequency, `f = ω / (2π)`.
///     * ``"k_complex"`` (:class:`complex`) — converged physical
///       wavenumber, `k = ω · √(μ₀ ε₀ ε(ω))`.
///     * ``"iterations"`` (int) — outer Newton iterations consumed.
///       Currently reported as ``max_iter`` on convergence and
///       ``max_iter`` on non-convergence; see Notes.
///     * ``"converged"`` (bool) — ``True`` iff the Newton fixed-point
///       relative-step test `|Δω/ω| < tol` triggered before
///       `max_iter`.
///
/// Notes
/// -----
/// `max_iter` and `tol` bound the **outer Newton fixed-point**
/// iteration only — the inner [`yee_fem::ComplexInverseIterEigen`] uses
/// its own defaults (1000 iter / 1e-8 Rayleigh-quotient tol), which are
/// orders of magnitude beyond what a user-facing Newton-iteration cap
/// would imply. The fine-grained inner-solver knobs remain available via
/// the Rust API (`DispersiveSolver::with_tuning`).
///
/// To report the precise Newton iteration count on convergence the binding
/// runs a small ladder: it retries
/// [`DispersiveSolver::solve_with_newton`] with `newton_max_iter` set to
/// successive integers from `1` up to `max_iter`, returning the first
/// successful result. This is O(`max_iter²`) inner solves on the
/// successful branch but `max_iter ≤ 8` keeps the overhead bounded. A
/// future revision can drop this ladder once `solve_with_newton` surfaces
/// the iteration counter on the success path natively (deferred per
/// ADR-0039 §"Material relocation").
///
/// The materials list is consumed in order — duplicate `tag` entries
/// register multiple materials under the same tag; the first match is
/// authoritative on lookup (see [`yee_fem::MaterialDatabase`]).
///
/// Raises
/// ------
/// ValueError
///     If any cavity extent is non-positive / non-finite, any subdivision
///     is zero, or any `materials` entry is malformed.
/// RuntimeError
///     If the inner sparse LU or the outer Newton fixed-point fails to
///     converge within `max_iter`.
#[pyfunction]
#[pyo3(signature = (a, b, d, nx, ny, nz, materials, omega_warm_start_hz, max_iter = 8, tol = 1e-6))]
#[allow(clippy::too_many_arguments)]
fn solve_cavity_dispersive<'py>(
    py: Python<'py>,
    a: f64,
    b: f64,
    d: f64,
    nx: usize,
    ny: usize,
    nz: usize,
    materials: &Bound<'py, PyList>,
    omega_warm_start_hz: f64,
    max_iter: usize,
    tol: f64,
) -> PyResult<Bound<'py, PyDict>> {
    // ---- 1. Build the cavity mesh -----------------------------------
    let mesh =
        TetMesh3D::cavity_uniform(a, b, d, nx, ny, nz).map_err(crate::errors::yee_mesh_to_py)?;

    // ---- 2. Parse the materials list into a MaterialDatabase --------
    let mut db = MaterialDatabase::new();
    for (i, item) in materials.iter().enumerate() {
        let mat_dict = item
            .cast::<PyDict>()
            .map_err(|e| PyValueError::new_err(format!("materials[{i}] must be a dict: {e}")))?;
        let (tag, mat) = material_from_pydict(mat_dict)?;
        db = db.with_material(tag, mat);
    }

    // ---- 3. Convert warm-start Hz → complex angular frequency -------
    let omega_warm_start = Complex64::new(2.0 * PI * omega_warm_start_hz, 0.0);

    // ---- 4. Configure the DispersiveSolver --------------------------
    //
    // Inner-solver tuning mirrors the D5 `dispersive_newton.rs` gate
    // tests: `max_iter = 1000` (the `DispersiveSolver::new` default) but
    // `tol = 1e-7` — one decade looser than the new-default 1e-8 because
    // the 10th deflated mode returned by `ComplexInverseIterEigen` sits
    // right at the working-precision boundary at 1e-8 for the WR-90
    // mesh density used here and in `dispersive_newton.rs`. The Python
    // `tol` argument bounds the *outer* Newton fixed-point only.
    let solver = DispersiveSolver::with_tuning(db, 1000, 1e-7);
    let max_iter = max_iter.max(1);

    // ---- 5. Solve --------------------------------------------------
    //
    // `sigma_factor = 2.5` mirrors the D4 / D5 fixture convention from
    // `dispersive_newton.rs` and the v0 `solve_cavity` shift heuristic.
    //
    // The ladder loop below retries `solve_with_newton` with successive
    // `newton_max_iter` values from `1` to `max_iter` so the precise
    // iteration count on a converged run can be reported back. On
    // non-convergence at the final attempt we fall through to the error
    // arm and surface `converged = false` with the last iterate.
    let mut omega_complex = Complex64::new(f64::NAN, f64::NAN);
    let mut k_complex = Complex64::new(f64::NAN, f64::NAN);
    let mut converged = false;
    let mut iterations = max_iter;
    for n in 1..=max_iter {
        let mut s = solver.clone();
        s.newton_max_iter = n;
        s.newton_tol = tol;
        match s.solve_with_newton(&mesh, omega_warm_start, 2.5) {
            Ok(eig) => {
                omega_complex = eig.omega;
                k_complex = eig.k_complex;
                converged = true;
                iterations = n;
                break;
            }
            Err(DispersiveError::NewtonDidNotConverge {
                last_omega,
                last_k_sq,
                last_residual: _,
            }) => {
                if n == max_iter {
                    // Final attempt failed — surface the last iterate
                    // so callers can diagnose.
                    omega_complex = last_omega;
                    k_complex = last_k_sq.sqrt();
                    converged = false;
                    iterations = max_iter;
                }
                // Otherwise keep climbing the ladder.
            }
            Err(DispersiveError::Underlying(e)) => return Err(crate::errors::yee_to_py(e)),
        }
    }

    // ---- 6. Pack into a Python dict --------------------------------
    let two_pi = 2.0 * PI;
    let f_complex = omega_complex / Complex64::new(two_pi, 0.0);

    let out = PyDict::new(py);
    out.set_item(
        "frequency_hz",
        PyComplex::from_doubles(py, f_complex.re, f_complex.im),
    )?;
    out.set_item(
        "k_complex",
        PyComplex::from_doubles(py, k_complex.re, k_complex.im),
    )?;
    out.set_item("iterations", iterations)?;
    out.set_item("converged", converged)?;
    Ok(out)
}
