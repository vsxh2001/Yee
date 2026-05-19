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

use nalgebra::Vector3;
use num_complex::Complex64;
use numpy::{IntoPyArray, PyArray2, PyArray3};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyComplex, PyDict, PyList};
use std::f64::consts::PI;
use yee_core::units::C0;
use yee_fem::{
    AbcOrder, DispersiveError, DispersiveSolver, FaceKind, FemEigenAssembly, InverseIterEigen,
    Material, MaterialDatabase, MaterialPole, OpenBoundarySolver, PortDefinition, SparseEigen,
};
use yee_mesh::{MaterialTag, TetMesh3D};

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(solve_cavity, m)?)?;
    m.add_function(wrap_pyfunction!(solve_cavity_dispersive, m)?)?;
    m.add_function(wrap_pyfunction!(solve_open_cavity, m)?)?;
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

/// Axis index from a Python axis string ("x" → 0, "y" → 1, "z" → 2).
fn axis_index(s: &str) -> PyResult<usize> {
    match s.to_ascii_lowercase().as_str() {
        "x" => Ok(0),
        "y" => Ok(1),
        "z" => Ok(2),
        other => Err(PyValueError::new_err(format!(
            "axis must be one of 'x', 'y', 'z' (got '{other}')"
        ))),
    }
}

/// Side enum from a Python side string ("low" → false, "high" → true).
fn side_is_high(s: &str) -> PyResult<bool> {
    match s.to_ascii_lowercase().as_str() {
        "low" => Ok(false),
        "high" => Ok(true),
        other => Err(PyValueError::new_err(format!(
            "side must be 'low' or 'high' (got '{other}')"
        ))),
    }
}

/// Per-port-face modal-profile descriptor.
///
/// Phase 4.fem.eig.2 v0 accepts a constant tangential E-field
/// ([`Self::Constant`]). Phase 4.fem.eig.3 F7 adds [`Self::Callable`],
/// wrapping a Python callable that takes a 3-tuple ``(x, y, z)`` and
/// returns a 3-tuple ``(Ex, Ey, Ez)`` so the wave-port can carry a
/// spatially-varying analytic profile (e.g. TE_{10} `ŷ · sin(π · x / a)`).
enum ModalProfile {
    /// Spatially-constant tangential `E_t` across every face centroid /
    /// Gauss point. The original v0 path.
    Constant(Vector3<f64>),
    /// Python callable `fn((x, y, z)) -> (Ex, Ey, Ez)` evaluated at each
    /// face centroid / Gauss point. `Py<PyAny>` is `Send + Sync` since
    /// PyO3 0.21; the call site re-acquires the GIL inside the closure
    /// via [`Python::attach`].
    Callable(Py<PyAny>),
}

impl ModalProfile {
    /// Cheap clone — duplicates the constant variant by `Copy` and the
    /// callable variant by [`Py::clone_ref`], which only increments the
    /// Python reference count.
    fn clone_with(&self, py: Python<'_>) -> Self {
        match self {
            Self::Constant(v) => Self::Constant(*v),
            Self::Callable(callable) => Self::Callable(callable.clone_ref(py)),
        }
    }

    /// Evaluate the modal profile at a world-space point.
    ///
    /// For the [`Self::Constant`] variant this returns the stored value
    /// unconditionally. For the [`Self::Callable`] variant this
    /// re-acquires the GIL, builds the input 3-tuple, calls the Python
    /// callable, and extracts the 3-tuple result. A panic propagates
    /// the error message via `expect`; callers must ensure the Python
    /// callable does not raise — the binding pre-validates this with a
    /// dry-run call in [`face_spec_from_pydict`].
    fn evaluate(&self, p: Vector3<f64>) -> Vector3<f64> {
        match self {
            Self::Constant(v) => *v,
            Self::Callable(callable) => Python::attach(|py| -> Vector3<f64> {
                let args = ((p.x, p.y, p.z),);
                let result = callable
                    .bind(py)
                    .call1(args)
                    .expect("modal_e_t callable raised during FEM assembly");
                let tup: (f64, f64, f64) = result
                    .extract()
                    .expect("modal_e_t callable returned non-(float, float, float) result");
                Vector3::new(tup.0, tup.1, tup.2)
            }),
        }
    }
}

/// Parsed face-plane descriptor produced by [`face_spec_from_pydict`].
///
/// Wave-port faces fill `port_id` and `modal_e_t`; ABC faces leave both
/// as `None`. `axis` is 0/1/2 for x/y/z and `side_high` is `false` for
/// the low coordinate of the cavity bounding box, `true` for the high.
struct FaceSpec {
    axis: usize,
    side_high: bool,
    port_id: Option<usize>,
    modal_e_t: Option<ModalProfile>,
}

/// Extract a [`FaceSpec`] from a Python face dict.
///
/// Shape:
/// `{"axis": "x"|"y"|"z", "side": "low"|"high", "port_id": int?,
///   "modal_e_t": (Ex, Ey, Ez) | Callable((x, y, z)) -> (Ex, Ey, Ez)?}`.
///
/// The callable form is Phase 4.fem.eig.3 F7 — required for the WR-90
/// thru-line gate which carries the analytic TE_{10} profile
/// `ŷ · sqrt(2/(a·b)) · sin(π·x/a)`. A constant 3-tuple still works as
/// the v0 fallback. The dispatcher tries `extract::<(f64, f64, f64)>`
/// first (cheap, no GIL re-acquisition at call time); on failure it
/// falls back to storing a `Py<PyAny>` that the closure re-binds via
/// [`Python::with_gil`].
fn face_spec_from_pydict(face_dict: &Bound<'_, PyDict>, want_port: bool) -> PyResult<FaceSpec> {
    let axis_str: String = face_dict
        .get_item("axis")?
        .ok_or_else(|| PyValueError::new_err("face dict missing required key 'axis'"))?
        .extract()
        .map_err(|e| PyValueError::new_err(format!("face 'axis' must be str: {e}")))?;
    let side_str: String = face_dict
        .get_item("side")?
        .ok_or_else(|| PyValueError::new_err("face dict missing required key 'side'"))?
        .extract()
        .map_err(|e| PyValueError::new_err(format!("face 'side' must be str: {e}")))?;
    let axis = axis_index(&axis_str)?;
    let side_high = side_is_high(&side_str)?;

    let (port_id, modal_e_t) = if want_port {
        let pid: usize = face_dict
            .get_item("port_id")?
            .ok_or_else(|| PyValueError::new_err("port face dict missing required key 'port_id'"))?
            .extract()
            .map_err(|e| PyValueError::new_err(format!("port_face 'port_id' must be int: {e}")))?;
        let modal_any = face_dict.get_item("modal_e_t")?.ok_or_else(|| {
            PyValueError::new_err("port face dict missing required key 'modal_e_t'")
        })?;
        let profile = if let Ok(tup) = modal_any.extract::<(f64, f64, f64)>() {
            ModalProfile::Constant(Vector3::new(tup.0, tup.1, tup.2))
        } else if modal_any.is_callable() {
            // Dry-run the callable on (0, 0, 0) so a malformed signature
            // surfaces a ValueError up-front rather than panicking inside
            // the FEM assembly. The Phase 4.fem.eig.3 F2 + F5 hot paths
            // call the profile many times per face; bailing on the first
            // call with `unwrap` keeps the closure trait-bounds clean
            // (the FEM closure type is `Fn(...) -> Vector3<f64>`, not
            // `Fn(...) -> Result<...>`).
            let probe = modal_any.call1(((0.0_f64, 0.0_f64, 0.0_f64),))?;
            let _: (f64, f64, f64) = probe.extract().map_err(|e| {
                PyValueError::new_err(format!(
                    "port_face 'modal_e_t' callable must return a (float, float, float) \
                     tuple: {e}"
                ))
            })?;
            ModalProfile::Callable(modal_any.unbind())
        } else {
            return Err(PyValueError::new_err(
                "port_face 'modal_e_t' must be either a tuple of 3 floats \
                 (Ex, Ey, Ez) or a callable taking a 3-tuple (x, y, z) and \
                 returning a 3-tuple (Ex, Ey, Ez)",
            ));
        };
        (Some(pid), Some(profile))
    } else {
        (None, None)
    };

    Ok(FaceSpec {
        axis,
        side_high,
        port_id,
        modal_e_t,
    })
}

/// Solve a swept open-boundary FEM driven analysis on a rectangular
/// cavity with ABC and modal wave-port faces.
///
/// Phase 4.fem.eig.2 step E6 + Phase 4.fem.eig.3 step F7 — Python
/// binding for [`yee_fem::OpenBoundarySolver::sweep`] (default) and
/// [`yee_fem::OpenBoundarySolver::sweep_matrix`] (with
/// ``multi_port=True``). Builds a Kuhn 6-tet uniform mesh via
/// [`yee_mesh::TetMesh3D::cavity_uniform`], constructs a
/// [`yee_fem::MaterialDatabase`] from the `materials` list (mirroring
/// :func:`yee.fem.solve_cavity_dispersive`), classifies every exterior
/// face by its centroid against the caller-supplied `port_faces` /
/// `abc_faces` axis-side specs (any face not on a tagged plane defaults
/// to PEC), assembles a constant-`β_mode = sqrt(k₀² − (π/a)²)` /
/// constant-`modal_e_t` wave-port descriptor per port, and runs the
/// per-frequency sparse-LU sweep.
///
/// Parameters
/// ----------
/// a, b, d : float
///     Cavity extents along the `x`, `y`, `z` axes, in metres.
/// nx, ny, nz : int
///     Brick subdivisions along each axis (Kuhn 6-tet decomposition).
/// materials : list[dict]
///     Same shape as :func:`yee.fem.solve_cavity_dispersive` — each
///     entry carries `tag`, `eps_inf`, `mu_r`, optional `poles`.
/// port_faces : list[dict]
///     Per-port-face descriptors of the shape
///     ``{"axis": "x"|"y"|"z", "side": "low"|"high", "port_id": int,
///        "modal_e_t": (Ex, Ey, Ez)}``. ``port_id`` indexes into the
///     wave-port descriptor list; multiple faces may share a port id.
///     ``modal_e_t`` is the **constant** tangential modal E-field at
///     every face centroid (Phase 4.fem.eig.2 v0 — per-Gauss-point
///     sampling is a 4.fem.eig.2.0.1 refinement).
/// abc_faces : list[dict]
///     ABC face descriptors of the shape
///     ``{"axis": "x"|"y"|"z", "side": "low"|"high"}``.
/// omegas_hz : list[float]
///     Real-valued sweep frequencies in Hz. Internally converted to
///     `ω = 2π · f` (rad/s).
/// coupled_whitney : bool, optional
///     If ``True``, the underlying solver wires the coupled
///     exact-Whitney-1 modal RHS + projection at three Gauss points
///     per port face (Phase 4.fem.eig.3 F1+F2). Default ``False``
///     reproduces the Phase 4.fem.eig.2 v2 + CCCCCCCCC lumped-centroid
///     behaviour bit-for-bit. Forwards to
///     [`yee_fem::OpenBoundarySolver::with_coupled_whitney`].
/// abc_order : str, optional
///     Order of the Engquist-Majda ABC bilinear form on
///     ``FaceKind::Abc``-tagged faces. One of ``"first"`` (default,
///     Phase 4.fem.eig.2 v2 1st-order Mur reproduced bit-for-bit) or
///     ``"second"`` (Phase 4.fem.eig.3 F4 2nd-order tangential-curl
///     correction, lowering the normal-incidence reflection floor from
///     ~-40 dB to ~-60 dB). Forwards to
///     [`yee_fem::OpenBoundarySolver::with_abc_order`]; any other
///     string raises ``ValueError``.
/// multi_port : bool, optional
///     If ``True``, the binding calls
///     [`yee_fem::OpenBoundarySolver::sweep_matrix`] (Phase 4.fem.eig.3
///     F5), returning the full `n_ports × n_ports` complex S-matrix
///     per frequency. Off-diagonal `S_{q,p}` entries are populated by
///     the per-excited-port column extraction of Sheen et al. 1990.
///     Default ``False`` falls back to the single-port
///     [`yee_fem::OpenBoundarySolver::sweep`] path; the returned
///     tensor still has shape `(n_omegas, n_ports, n_ports)` with
///     off-diagonals set to zero.
///
/// Returns
/// -------
/// numpy.ndarray
///     `(n_omegas, n_ports, n_ports)` complex128 S-parameter tensor.
///     With ``multi_port=False`` (default): single-port
///     S-parameters only; off-diagonal `S_{p,q}` entries for `p ≠ q`
///     are zero. With ``multi_port=True``: full multi-port S-matrix
///     per frequency.
///
/// Notes
/// -----
/// The `β_mode(ω) = sqrt((ω/c)² − (π/a)²)` closure is the rectangular
/// TE_{10} dispersion against the **first** cavity extent `a`. Callers
/// driving a non-rectangular cross-section should evaluate the modal
/// dispersion outside this binding and call the underlying Rust API
/// directly.
///
/// Below-cutoff frequencies (`ω/c < π/a`) yield `β_mode = 0`; the
/// per-frequency assembly still completes, but the wave-port modal
/// source vanishes there.
///
/// Raises
/// ------
/// ValueError
///     On any malformed cavity extent / subdivision / materials entry
///     / face descriptor, if any tagged axis-side plane is not an
///     actual exterior face of the mesh, or if ``abc_order`` is not
///     one of ``"first"`` / ``"second"``.
/// RuntimeError
///     If any per-frequency sparse LU fails to factor.
#[pyfunction]
#[pyo3(signature = (
    a, b, d, nx, ny, nz, materials, port_faces, abc_faces, omegas_hz,
    coupled_whitney = false, abc_order = "first", multi_port = false
))]
#[allow(clippy::too_many_arguments)]
fn solve_open_cavity<'py>(
    py: Python<'py>,
    a: f64,
    b: f64,
    d: f64,
    nx: usize,
    ny: usize,
    nz: usize,
    materials: &Bound<'py, PyList>,
    port_faces: &Bound<'py, PyList>,
    abc_faces: &Bound<'py, PyList>,
    omegas_hz: &Bound<'py, PyList>,
    coupled_whitney: bool,
    abc_order: &str,
    multi_port: bool,
) -> PyResult<Bound<'py, PyArray3<Complex64>>> {
    // ---- 1. Build the cavity mesh -----------------------------------
    let mesh =
        TetMesh3D::cavity_uniform(a, b, d, nx, ny, nz).map_err(crate::errors::yee_mesh_to_py)?;

    // ---- 2. Parse materials → MaterialDatabase ----------------------
    let mut db = MaterialDatabase::new();
    for (i, item) in materials.iter().enumerate() {
        let mat_dict = item
            .cast::<PyDict>()
            .map_err(|e| PyValueError::new_err(format!("materials[{i}] must be a dict: {e}")))?;
        let (tag, mat) = material_from_pydict(mat_dict)?;
        db = db.with_material(tag, mat);
    }

    // ---- 3. Parse abc_faces + port_faces specs ----------------------
    //
    // Each entry resolves to a `(axis, side_high)` tuple identifying a
    // half-space plane of the cavity bounding box. Wave-port entries
    // also carry a port_id (into the `ports` vector below) and a
    // constant modal_e_t vector.
    let mut abc_planes: Vec<(usize, bool)> = Vec::with_capacity(abc_faces.len());
    for (i, item) in abc_faces.iter().enumerate() {
        let face_dict = item
            .cast::<PyDict>()
            .map_err(|e| PyValueError::new_err(format!("abc_faces[{i}] must be a dict: {e}")))?;
        let spec = face_spec_from_pydict(face_dict, false)?;
        abc_planes.push((spec.axis, spec.side_high));
    }

    // Wave-port plane list — each entry carries axis/side/port_id and a
    // `ModalProfile` (constant tuple or Python callable). We also
    // collect the unique port_id set and build one `PortDefinition` per
    // id with that port's modal profile (taking the first face's
    // profile per id — v0 single-modal-source-per-port convention).
    let mut port_planes: Vec<(usize, bool, usize)> = Vec::with_capacity(port_faces.len());
    let mut max_port_id: i64 = -1;
    let mut port_modal_e_t: std::collections::HashMap<usize, ModalProfile> =
        std::collections::HashMap::new();
    for (i, item) in port_faces.iter().enumerate() {
        let face_dict = item
            .cast::<PyDict>()
            .map_err(|e| PyValueError::new_err(format!("port_faces[{i}] must be a dict: {e}")))?;
        let spec = face_spec_from_pydict(face_dict, true)?;
        let port_id = spec
            .port_id
            .expect("face_spec_from_pydict want_port=true returns Some");
        let modal_profile = spec
            .modal_e_t
            .expect("face_spec_from_pydict want_port=true returns Some");
        port_planes.push((spec.axis, spec.side_high, port_id));
        if (port_id as i64) > max_port_id {
            max_port_id = port_id as i64;
        }
        port_modal_e_t.entry(port_id).or_insert(modal_profile);
    }
    let n_ports = if max_port_id < 0 {
        0
    } else {
        (max_port_id + 1) as usize
    };

    // ---- 4. Build PortDefinitions -----------------------------------
    //
    // `β_mode(ω) = sqrt((ω/c)² − (π/a)²)` (TE_{10} dispersion on the
    // first cavity extent). Below-cutoff returns 0 — same convention as
    // `yee_validation::run_fem_eig_003_wr90_stub_abc`.
    //
    // `modal_e_t(x) = profile.evaluate(x)` — Phase 4.fem.eig.3 F7
    // extends the v0 constant-profile path so the closure dispatches to
    // either a stored `Vector3<f64>` (constant) or a Python callable
    // (re-acquires the GIL inside the closure). The
    // [`ModalProfile::Callable`] path is what enables the WR-90
    // thru-line gate to carry the analytic TE_{10} profile `ŷ · sin(π
    // · x / a)`.
    let mut ports: Vec<PortDefinition> = Vec::with_capacity(n_ports);
    let a_for_beta = a;
    for p in 0..n_ports {
        let profile = port_modal_e_t
            .get(&p)
            .map(|prof| prof.clone_with(py))
            .ok_or_else(|| {
                PyValueError::new_err(format!(
                    "port_faces does not register any face for port_id = {p} \
                     (port ids must be a contiguous range [0, n_ports))"
                ))
            })?;
        ports.push(PortDefinition {
            beta_mode: Box::new(move |omega: f64| -> f64 {
                let k0_sq = (omega / C0).powi(2);
                let kc_sq = (PI / a_for_beta).powi(2);
                let arg = k0_sq - kc_sq;
                if arg <= 0.0 { 0.0 } else { arg.sqrt() }
            }),
            modal_e_t: Box::new(move |p: Vector3<f64>| -> Vector3<f64> { profile.evaluate(p) }),
        });
    }

    // ---- 5. Classify exterior faces by centroid ---------------------
    //
    // Build a placeholder all-PEC solver to recover the canonical
    // exterior-face centroid ordering (mirror of the
    // `yee_validation::run_fem_eig_003_wr90_stub_abc` pattern). Then
    // tag each centroid against the abc/port plane lists.
    //
    // Count exterior faces by walking the tet face-incidence map (one
    // tet face per `TET_FACES` entry; faces with multiplicity 1 are
    // exterior). This mirrors the helper in the yee-validation driver.
    let n_exterior = {
        let mut face_map: std::collections::HashMap<[usize; 3], usize> =
            std::collections::HashMap::new();
        const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        for tet in &mesh.tetrahedra {
            for &[ai, bi, ci] in TET_FACES.iter() {
                let mut key = [tet[ai], tet[bi], tet[ci]];
                key.sort_unstable();
                *face_map.entry(key).or_insert(0) += 1;
            }
        }
        face_map.values().filter(|&&c| c == 1).count()
    };
    let placeholder = OpenBoundarySolver::new(
        &mesh,
        vec![FaceKind::Pec; n_exterior],
        Vec::new(),
        MaterialDatabase::new(),
    )
    .map_err(crate::errors::yee_to_py)?;
    let centroids = placeholder.exterior_face_centroids();

    // Plane positions: `low` ⇒ coordinate 0, `high` ⇒ coordinate
    // (a, b, d)[axis]. The cavity_uniform mesh spans `[0, a] × [0, b] ×
    // [0, d]` exactly.
    let extents = [a, b, d];
    let plane_tol = 1e-9;
    let plane_value =
        |axis: usize, side_high: bool| -> f64 { if side_high { extents[axis] } else { 0.0 } };

    let mut face_kinds: Vec<FaceKind> = Vec::with_capacity(centroids.len());
    for c in &centroids {
        let coord = [c.x, c.y, c.z];
        // Wave-port classification first — takes precedence over PEC
        // default, but the underlying solver enforces PEC-precedence on
        // shared edges if a wave-port face neighbours a PEC face.
        let mut kind = FaceKind::Pec;
        for &(axis, side_high, port_id) in &port_planes {
            let target = plane_value(axis, side_high);
            if (coord[axis] - target).abs() < plane_tol {
                kind = FaceKind::WavePort(port_id);
                break;
            }
        }
        if matches!(kind, FaceKind::Pec) {
            for &(axis, side_high) in &abc_planes {
                let target = plane_value(axis, side_high);
                if (coord[axis] - target).abs() < plane_tol {
                    kind = FaceKind::Abc;
                    break;
                }
            }
        }
        face_kinds.push(kind);
    }

    // Sanity check: every tagged plane must intersect at least one
    // exterior face. Empty intersection usually means the caller passed
    // the wrong axis / side string.
    for &(axis, side_high, _) in &port_planes {
        let target = plane_value(axis, side_high);
        let any = centroids
            .iter()
            .any(|c| (([c.x, c.y, c.z])[axis] - target).abs() < plane_tol);
        if !any {
            return Err(PyValueError::new_err(format!(
                "solve_open_cavity: no exterior face on port plane axis={axis} side={} \
                 (target coord = {target})",
                if side_high { "high" } else { "low" }
            )));
        }
    }
    for &(axis, side_high) in &abc_planes {
        let target = plane_value(axis, side_high);
        let any = centroids
            .iter()
            .any(|c| (([c.x, c.y, c.z])[axis] - target).abs() < plane_tol);
        if !any {
            return Err(PyValueError::new_err(format!(
                "solve_open_cavity: no exterior face on abc plane axis={axis} side={} \
                 (target coord = {target})",
                if side_high { "high" } else { "low" }
            )));
        }
    }

    // ---- 6. Parse the abc_order kwarg → AbcOrder enum ---------------
    //
    // Case-insensitive match; any other string surfaces as ValueError
    // per the public docstring. Mirrors the axis / side parser pattern
    // above.
    let abc_order_enum = match abc_order.to_ascii_lowercase().as_str() {
        "first" => AbcOrder::First,
        "second" => AbcOrder::Second,
        other => {
            return Err(PyValueError::new_err(format!(
                "solve_open_cavity: abc_order must be 'first' or 'second' \
                 (got '{other}')"
            )));
        }
    };

    // ---- 7. Build the real solver + run the sweep -------------------
    //
    // Phase 4.fem.eig.3 F7: `coupled_whitney` and `abc_order` are
    // applied via the F2 / F4 builder methods on `OpenBoundarySolver`.
    // Defaults (`false` / `AbcOrder::First`) reproduce the v2 +
    // CCCCCCCCC behaviour bit-for-bit.
    let solver = OpenBoundarySolver::new(&mesh, face_kinds, ports, db)
        .map_err(crate::errors::yee_to_py)?
        .with_coupled_whitney(coupled_whitney)
        .with_abc_order(abc_order_enum);

    let omegas: Vec<f64> = {
        let mut out: Vec<f64> = Vec::with_capacity(omegas_hz.len());
        for (i, item) in omegas_hz.iter().enumerate() {
            let f: f64 = item.extract().map_err(|e| {
                PyValueError::new_err(format!("omegas_hz[{i}] must be a float: {e}"))
            })?;
            out.push(2.0 * PI * f);
        }
        out
    };
    if omegas.is_empty() {
        return Err(PyValueError::new_err(
            "solve_open_cavity: omegas_hz must contain at least one frequency",
        ));
    }
    if n_ports == 0 {
        return Err(PyValueError::new_err(
            "solve_open_cavity: at least one wave-port face is required \
             (no driver → no S-parameters to compute)",
        ));
    }

    // ---- 8. Pack into a (n_omegas, n_ports, n_ports) tensor ---------
    //
    // `multi_port = False` → single-port `sweep` path; only the
    // diagonal `S_{p,p}` entries are populated and the off-diagonals
    // default to zero (Phase 4.fem.eig.2 v0 shape contract).
    //
    // `multi_port = True` → Phase 4.fem.eig.3 F5 `sweep_matrix` path;
    // every entry of the per-frequency `n_ports × n_ports` matrix is
    // computed via per-excited-port column extraction.
    let n_omegas = omegas.len();
    let mut buf: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n_omegas * n_ports * n_ports];
    if multi_port {
        let sweep = solver
            .sweep_matrix(&omegas)
            .map_err(crate::errors::yee_to_py)?;
        for (k, s_k) in sweep.s.iter().enumerate().take(n_omegas) {
            for q in 0..n_ports {
                for p in 0..n_ports {
                    let idx = (k * n_ports + q) * n_ports + p;
                    buf[idx] = s_k[(q, p)];
                }
            }
        }
    } else {
        let sweep = solver.sweep(&omegas).map_err(crate::errors::yee_to_py)?;
        for (p, port_sweep) in sweep.s_pp.iter().enumerate().take(n_ports) {
            for (k, &s_pp_k) in port_sweep.iter().enumerate().take(n_omegas) {
                let idx = (k * n_ports + p) * n_ports + p;
                buf[idx] = s_pp_k;
            }
        }
    }
    let arr = numpy::ndarray::Array3::from_shape_vec((n_omegas, n_ports, n_ports), buf)
        .expect("buffer length matches (n_omegas, n_ports, n_ports) by construction");
    Ok(arr.into_pyarray(py))
}
