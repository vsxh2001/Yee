//! Source helpers for the FDTD walking skeleton.
//!
//! Phase 2.0 shipped a single point-source primitive: a Gaussian-in-time pulse
//! added (soft source) to a chosen cell of `E_z`. Phase 2.fdtd.5 adds a
//! total-field / scattered-field (TF/SF) plane-wave source, see
//! [`PlaneWaveSource`].
//!
//! Hard sources, modal sources, and lumped ports remain Phase 2.1+ work.
//!
//! ## Phase 2.fdtd.5.2 design notes — j/k-face SF corrections
//!
//! **Assumption being challenged:** Phase 2.fdtd.5.1 shipped only `i`-face
//! TF/SF corrections, on the reasoning that for a `+x` `E_z`-polarized
//! plane wave the incident `H_inc_x`, `H_inc_z`, `E_inc_x`, `E_inc_y` are
//! all zero, so the j- and k-face stencils carry no spurious **incident**
//! contribution. That argument is correct for the *incident* leg but
//! misses the *scattered* leg: the j-face has a TF-vs-SF discontinuity
//! in `E_z` (it equals `E_inc_z` plus scattered inside the box, just
//! scattered outside), and the standard Yee `H_x` update at `j = j0 - 1`
//! and `j = j1` straddles that discontinuity in its `∂E_z / ∂y` term,
//! mixing TF and SF `E_z` and emitting spurious scattered field into the
//! SF region. Symmetrically, the `E_x` update at `k = k0` and
//! `k = k1 + 1` straddles the z-discontinuity in `H_y` via its
//! `∂H_y / ∂z` term. With those four corrections added, the finite-box
//! configuration's TF/SF contrast jumps from ~6× (Phase 2.fdtd.5.1
//! empirical pin) to >100× (Phase 2.fdtd.5.2 target), and slab
//! geometry — which puts the j/k faces inside CPML — is unaffected.
//!
//! **Which curl stencils need a correction for `+x` `E_z` polarization**
//! (only `E_inc_z(x)` and `H_inc_y(x)` are non-zero):
//!
//! - `H_y` curl has `∂E_z / ∂x` — i-face straddle. (5.1, shipped.)
//! - `H_x` curl has `∂E_z / ∂y` — j-face straddle. (5.2, this commit.)
//! - `E_z` curl has `∂H_y / ∂x` — i-face straddle. (5.1, shipped.)
//! - `E_x` curl has `∂H_y / ∂z` — k-face straddle. (5.2, this commit.)
//!
//! All other E/H components' curls involve only zero-incident pairs
//! (`E_inc_x = E_inc_y = H_inc_x = H_inc_z = 0`), so they need no
//! correction for `+x` `E_z` polarization. Arbitrary-polarization /
//! oblique-incidence support lands in Phase 2.fdtd.5.3+.
//!
//! **Approach:** extend the existing i-face apply pattern to j and k
//! faces. Each face/component pair is one inclusive 2-D loop matching
//! the i-face's `j0..=j1, k0..k_hi` index conventions. The sign of the
//! correction comes from the side of the box and the orientation of the
//! straddled discontinuity, derived in the per-face comment blocks
//! below.
//!
//! **Reference:** Taflove & Hagness, *Computational Electrodynamics*
//! (3rd ed.) §5.10 — 3-D TF/SF for a rectangular Huygens surface.
//!
//! **DoD:** the `tests/plane_wave_finite_box.rs` contrast must rise
//! from ~6× to ≥ 100×, and the slab variant
//! (`tests/plane_wave_propagation.rs`) must not regress below its
//! previous ~2676×.
//!
//! ## Phase 2.fdtd.5.3 design notes — oblique incidence + general polarization
//!
//! **Assumption being challenged:** Phase 2.fdtd.5/5.1/5.2 hard-coded
//! `k_hat = +x_hat` with `E` along `z_hat` and `H` along `y_hat`. Three
//! of the six potentially-incident field components were therefore
//! identically zero, which collapsed the general 12-face-correction
//! Taflove §5.10 stencil down to four (i, j-faces for H; i, k-faces
//! for E). This track removes that restriction.
//!
//! ### 1-D auxiliary incident-field grid along `k_hat`
//!
//! For an arbitrary unit propagation vector
//! ```text
//!   k_hat = (sin θ cos φ, sin θ sin φ, cos θ)
//! ```
//! the incident field is `E_inc(r, t) = E_inc_vec · f(t − r·k_hat/c)`,
//! a plane wave whose phase varies only along `k_hat`. We therefore
//! keep a single **1-D Yee grid** along `k_hat` carrying the *scalar*
//! incident amplitude `f(s, t)` (where `s` is signed distance along
//! `k_hat`) and update it in lockstep with the 3-D solver using the
//! same 1-D leapfrog as the 5.x normal-incidence path. The 3-D vector
//! decomposition happens **at the box-face correction stencil**, not
//! in the aux grid.
//!
//! At each face Yee node `r` we evaluate
//! ```text
//!   s(r) = (r − r_ref) · k_hat
//! ```
//! where `r_ref` is the TF-box corner closest to the wave source
//! (for `k_hat` with non-negative components, this is the
//! `(i0·dx, j0·dy, k0·dz)` corner). `s ≥ 0` for every TF box node
//! when `k_hat` has non-negative components; the aux grid is sized
//! `n = ceil(S_max/dx) + 2·pad + 1` where `S_max` is the projected
//! box diagonal `Σ_α (i1_α − i0_α)·dx_α·|k_hat_α|`. Index `pad`
//! corresponds to `s = 0` (the TF reference corner); the source
//! cell is at index `0` (one `dx · pad` upstream), and the Mur ABC
//! sits at the far end.
//!
//! Sampling: linear interpolation in `s` into `inc_e` (or `inc_h`,
//! which lives at `s + dx/2`). At normal incidence this collapses
//! to the 5.2 lookup with zero interpolation error; off-axis the
//! interpolation introduces `O((dx · sin θ)²)` dispersion / phase
//! mismatch, which is acceptable inside the test tolerance.
//!
//! ### Polarization parametrization
//!
//! Inputs are three angles `(θ, φ, ψ)`. `(θ, φ)` define `k_hat` as
//! above; `ψ ∈ [0, 2π)` rotates the E-vector inside the plane
//! perpendicular to `k_hat`. Build the spherical basis
//! ```text
//!   e_theta_hat = ( cos θ cos φ,  cos θ sin φ, −sin θ )
//!   e_phi_hat   = (−sin φ,        cos φ,        0      )
//! ```
//! Then
//! ```text
//!   E_inc_hat = cos ψ · e_theta_hat + sin ψ · e_phi_hat
//!   H_inc_hat = −(k_hat × E_inc_hat)
//!             = −(cos ψ · e_phi_hat − sin ψ · e_theta_hat)
//! ```
//! (using `k̂ × e_θ̂ = e_φ̂`, `k̂ × e_φ̂ = −e_θ̂`.) The minus sign on
//! `H_inc_hat` absorbs the 1-D Yee aux grid's intrinsic sign
//! convention: the 1-D leapfrog gives `inc_h ≈ −inc_e/η₀` for a
//! `+s`-propagating wave (i.e. `inc_h` carries the negative of the
//! physical H magnitude). With `H_inc_hat = −(k̂×E_inc_hat)`, the
//! reconstruction
//! ```text
//!   E_inc_α(r, t) = E_inc_hat_α · inc_e[ interp(s(r)) ]
//!   H_inc_α(r, t) = H_inc_hat_α · inc_h[ interp(s(r)) ]
//! ```
//! correctly yields the Maxwell-consistent H field (in particular,
//! for k̂=+x̂, E along +ẑ — the 5.2 legacy case — this gives
//! H_inc_hat = +ŷ, and combined with the negative `inc_h` produces
//! the physical Maxwell `H_y = −E_z/η₀ < 0`).
//! Both `E_inc_hat` and `H_inc_hat` are dimensionless unit vectors;
//! `inc_e` and `inc_h` are scalar.
//!
//! ### Why the six-face structure stays the same
//!
//! The TF/SF correction stencils at each box face are derived purely
//! from "this curl term straddled the boundary; subtract the
//! incident contribution it picked up." That derivation does not
//! depend on the propagation direction — only on (a) which curl
//! terms have a non-zero incident component, and (b) what that
//! incident value is at the face node. For arbitrary polarization
//! **all twelve face-component pairs are potentially nonzero**:
//!
//! - i-faces: H_y reads E_inc_z, H_z reads E_inc_y, E_y reads H_inc_z,
//!   E_z reads H_inc_y.
//! - j-faces: H_x reads E_inc_z, H_z reads E_inc_x, E_x reads H_inc_z,
//!   E_z reads H_inc_x.
//! - k-faces: H_x reads E_inc_y, H_y reads E_inc_x, E_x reads H_inc_y,
//!   E_y reads H_inc_x.
//!
//! The 5.2 path was the special case where six of these twelve
//! incident references collapsed to zero (E_inc_x, E_inc_y, H_inc_x,
//! H_inc_z all identically zero). The 5.3 implementation populates
//! all twelve with `(E_inc_hat_α · inc_e[…])` / `(H_inc_hat_α · inc_h[…])`
//! from the aux grid + projection at face-node coordinates.
//!
//! ### Back-compat
//!
//! [`PlaneWaveSource::new`] continues to construct a normal-incidence
//! `(θ=0, φ=0, ψ=0)` source — the 5.2 path — bit-for-bit. The new
//! [`PlaneWaveSource::with_oblique_incidence`] takes the three angles
//! and drives the general kernel. The implementation specializes the
//! `θ=0, φ=0, ψ=0` case at construction (sets all polarization
//! components to their normal-incidence values) but uses the same
//! 12-face correction kernel; the normal-incidence regression test
//! pins the resulting contrast to within 1% of the 5.2 floor.
//!
//! **Reference:** Taflove & Hagness, *Computational
//! Electrodynamics* (3rd ed.) §5.10 (3-D TF/SF for arbitrary plane
//! wave on a rectangular Huygens surface).
//!
//! **DoD:** (a) the 5.2 finite-box contrast at `(θ=0, φ=0, ψ=0)`
//! must remain within 1% of its previous ≥ 1e14× floor; (b) at
//! `(θ=30°, φ=45°)` with `E_phi` polarization the finite-box contrast
//! must clear 1000× (the loose interpolation-limited gate); (c)
//! `θ = 85°` must run without NaN / panic, even if contrast degrades.
//!
//! ## Phase 2.fdtd.5.3.1 — dispersion-matched 1-D auxiliary step
//!
//! Phase 2.fdtd.5.3 shipped the 12-face kernel with `ds_aux = dx`. At
//! oblique 30°/45° this plateaus at ~14.5× contrast — well below the
//! 1000× DoD. The leakage is **dispersion mismatch**: the 1-D Yee aux
//! propagates at the 1-D phase velocity (a function of `ω·dt` and
//! `ds_aux`) while the 3-D solver propagates the same plane wave at
//! the *3-D* numerical phase velocity along `k_hat`. The two velocities
//! agree exactly on-axis (5.2 case, `k_hat = +x̂` and `ds_aux = dx`)
//! but disagree off-axis. Across an ~30-cell projected box diagonal
//! the accumulated phase mismatch drives the per-face leakage to
//! ~5-10%, dwarfing the per-cell linear interpolation error.
//!
//! Taflove §5.10.5's remedy: choose `ds_aux` such that the 1-D Yee
//! numerical wavenumber at the source frequency equals the projected
//! 3-D numerical wavenumber along `k_hat`. The 1-D Yee dispersion
//! relation is
//! ```text
//!   sin(k_1D · ds_aux / 2) / ds_aux = sin(ω·dt / 2) / (c·dt)
//! ```
//! and the 3-D Yee dispersion relation (cubic cells) is
//! ```text
//!   sin²(k · k̂_x · dx/2) + sin²(k · k̂_y · dx/2) + sin²(k · k̂_z · dx/2)
//!       = (dx/(c·dt))² · sin²(ω·dt/2)
//! ```
//! Solving the 3-D equation for `k = k_3D` (numerical wavenumber along
//! `k_hat`) and substituting `k_1D = k_3D` gives a transcendental
//! equation in `ds_aux`:
//! ```text
//!   f(ds) = sin(k_3D · ds / 2) / ds − sin(ω·dt/2) / (c·dt) = 0
//! ```
//! `f` has a degenerate root at `ds → 0` that an unguarded Newton
//! or fixed-point solver collapses to. The physical root sits within
//! a few percent of `dx`. **Bisection on `[0.5·dx, 2.0·dx]`** is the
//! correct discriminator: deterministic, immune to the trivial-root
//! attractor, converges in O(50) iterations to `|f| < 1e-12`.
//!
//! `compute_aux_step` (private) implements both bisection steps (3-D `k`
//! solve, then 1-D `ds_aux` solve). [`PlaneWaveSource::with_oblique_incidence`]
//! invokes it; the new
//! [`PlaneWaveSource::with_oblique_incidence_match`] takes an explicit
//! `dispersion_match: bool` flag and is used by the sanity test to
//! reproduce the Phase 2.fdtd.5.3 14.5× contrast for back-compat
//! comparison.
//!
//! **DoD:** the 30°/45° contrast must clear 1000× (the original
//! Phase 2.fdtd.5.3 DoD that was deferred). The normal-incidence
//! regression must stay at its 1e14× floor.

use std::f64::consts::{PI, TAU};

use yee_core::units::{C0, EPS0, MU0};

use crate::grid::YeeGrid;

/// Solve the 3-D Yee numerical dispersion relation for the
/// wavenumber magnitude `k` along the unit propagation vector
/// `k_hat = (sin θ cos φ, sin θ sin φ, cos θ)` at angular frequency
/// `ω`, given cubic cell size `dx` and time step `dt`.
///
/// The relation (cubic Yee cells) is
/// ```text
///   sin²(k · k̂_x · dx/2) + sin²(k · k̂_y · dx/2) + sin²(k · k̂_z · dx/2)
///       = (dx / (c·dt))² · sin²(ω·dt/2)
/// ```
/// Solved by bisection on `k ∈ (0, π/(max|k̂_α|·dx))`. Returns the
/// physical (smallest-positive) root.
fn solve_k_3d(theta: f64, phi: f64, dx: f64, c0: f64, omega: f64, dt: f64) -> f64 {
    let (sin_t, cos_t) = theta.sin_cos();
    let (sin_p, cos_p) = phi.sin_cos();
    let kh = [sin_t * cos_p, sin_t * sin_p, cos_t];

    // RHS of the dispersion equation: (dx/(c·dt))² · sin²(ω·dt/2).
    let s_omega = (omega * dt / 2.0).sin();
    let rhs = (dx / (c0 * dt)).powi(2) * s_omega * s_omega;

    // LHS as a function of k:  Σ sin²(k · k̂_α · dx/2).
    let lhs = |k: f64| -> f64 {
        let mut s = 0.0;
        for &kha in &kh {
            let a = (k * kha * dx / 2.0).sin();
            s += a * a;
        }
        s
    };

    // f(k) = LHS(k) - RHS. f(0) = -rhs < 0.  LHS is strictly increasing
    // in k on (0, π/(max|k̂_α|·dx)), so f crosses zero at most once.
    let f = |k: f64| -> f64 { lhs(k) - rhs };

    // Bracket: lower = small ε·ω/c, upper just below π/(max|k̂_α|·dx).
    // Use max_kh as the binding component; if all components are zero
    // (impossible for unit k_hat, but guard anyway) fall back to k = ω/c.
    let max_kh = kh.iter().fold(0.0_f64, |acc, &v| acc.max(v.abs()));
    if max_kh <= 0.0 {
        return omega / c0;
    }
    let mut lo = 1.0e-12_f64.max(omega / (10.0 * c0));
    let mut hi = (PI / (max_kh * dx)) * 0.999999;
    // Make sure f(lo) < 0 < f(hi). For physical inputs (sub-Nyquist
    // ω·dt) this holds, but tighten / widen as needed.
    let mut f_lo = f(lo);
    let mut f_hi = f(hi);
    // If f(lo) is already ≥ 0 (extremely small ω), shrink lo.
    let mut tries = 0;
    while f_lo >= 0.0 && tries < 16 {
        lo *= 0.1;
        f_lo = f(lo);
        tries += 1;
    }
    tries = 0;
    while f_hi <= 0.0 && tries < 16 {
        // Push hi just a touch closer to the singular limit. If we hit
        // it without bracketing, return the analytical guess as a
        // safe fallback.
        hi = lo + (hi - lo) * 0.5;
        f_hi = f(hi);
        tries += 1;
        if hi <= lo {
            break;
        }
    }
    if !(f_lo < 0.0 && f_hi > 0.0) {
        // Unbracketed — return the continuum wavenumber as a fallback.
        return omega / c0;
    }

    // Bisect to |f| < 1e-14 or 80 iterations.
    let _ = (f_lo, f_hi); // bracketing checks done above; no longer needed
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let f_mid = f(mid);
        if f_mid.abs() < 1.0e-14 {
            return mid;
        }
        if f_mid < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo) < 1.0e-18 {
            break;
        }
    }
    0.5 * (lo + hi)
}

/// Compute the dispersion-matched 1-D auxiliary-grid step `ds_aux` so
/// that the 1-D Yee leapfrog (driving the TF/SF incident field) propagates
/// at the same numerical phase velocity along `k_hat` as the 3-D Yee
/// grid does at angular frequency `ω`.
///
/// Implements Taflove & Hagness, *Computational Electrodynamics* (3rd
/// ed.) §5.10.5. The 1-D Yee dispersion relation is
/// ```text
///   sin(k_1D · ds_aux / 2) / ds_aux = sin(ω·dt / 2) / (c·dt)
/// ```
/// We choose `ds_aux` so `k_1D = k_3D` (the projected 3-D numerical
/// wavenumber along `k_hat`); substituting gives the transcendental
/// equation
/// ```text
///   f(ds) = sin(k_3D · ds / 2) / ds − sin(ω·dt/2) / (c·dt) = 0
/// ```
/// `f` has a degenerate root at `ds → 0` (because `sin(x)/x → 1` faster
/// than `x` shrinks, the limit of `sin(k·ds/2)/ds` is `k/2`, which is
/// generally not equal to the RHS — so the singularity is at `ds = 0`
/// not the root we want). The physical root sits within a few percent
/// of `dx`. We use **bisection** on `[0.5·dx, 2.0·dx]` (widened to
/// `[0.1·dx, 5·dx]` if the initial bracket fails to straddle zero),
/// which is deterministic, robust to the trivial-root attractor, and
/// converges in O(50) iterations to `|f| < 1e-12`.
///
/// On bracketing failure (extremely exotic geometry or ill-conditioned
/// inputs), falls back to `ds_aux = dx` and emits a `tracing::warn!`
/// rather than panicking.
fn compute_aux_step(theta: f64, phi: f64, dx: f64, c0: f64, omega: f64, dt: f64) -> f64 {
    let k_3d = solve_k_3d(theta, phi, dx, c0, omega, dt);

    let rhs = (omega * dt / 2.0).sin() / (c0 * dt);
    let f = |ds: f64| -> f64 { (k_3d * ds / 2.0).sin() / ds - rhs };

    // Try primary bracket [0.5·dx, 2.0·dx] then widen to [0.1·dx, 5·dx]
    // if the function doesn't straddle zero across the first.
    let brackets = [(0.5 * dx, 2.0 * dx), (0.1 * dx, 5.0 * dx)];
    for &(lo0, hi0) in &brackets {
        let mut lo = lo0;
        let mut hi = hi0;
        let mut f_lo = f(lo);
        let mut f_hi = f(hi);
        if f_lo == 0.0 {
            return lo;
        }
        if f_hi == 0.0 {
            return hi;
        }
        if f_lo.signum() == f_hi.signum() {
            continue; // try a wider bracket
        }

        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            let f_mid = f(mid);
            if f_mid.abs() < 1.0e-12 {
                return mid;
            }
            if f_mid.signum() == f_lo.signum() {
                lo = mid;
                f_lo = f_mid;
            } else {
                hi = mid;
                f_hi = f_mid;
            }
            if (hi - lo) < 1.0e-15 * dx {
                break;
            }
        }
        let root = 0.5 * (lo + hi);
        if root > 0.0 && root.is_finite() {
            return root;
        }
    }

    tracing::warn!(
        target: "yee_fdtd::sources",
        theta = theta,
        phi = phi,
        dx = dx,
        omega = omega,
        dt = dt,
        k_3d = k_3d,
        "compute_aux_step: bisection failed to bracket the dispersion-matched root; falling back to ds_aux = dx"
    );
    dx
}

/// 4-point cubic Lagrange interpolation on a uniformly-spaced 1-D array
/// (the 1-D auxiliary TF/SF incident-field grid). `f` is the
/// fractional index (with `arr[i]` at integer `i`). Falls back to
/// linear interpolation in the boundary cells (where the 4-point
/// stencil would step out of bounds), and clamps to the array end
/// values outside `[0, n-1]`.
///
/// The cubic-Lagrange basis evaluated at the four samples
/// `arr[m-1], arr[m], arr[m+1], arr[m+2]` (where `m = floor(f)`,
/// `t = f - m ∈ [0,1)`) is
/// ```text
///   L_{-1}(t) = -t(t-1)(t-2)/6
///   L_0(t)    =  (t+1)(t-1)(t-2)/2
///   L_1(t)    = -(t+1)t(t-2)/2
///   L_2(t)    =  (t+1)t(t-1)/6
/// ```
fn sample_aux_cubic(arr: &[f64], f: f64) -> f64 {
    let n = arr.len();
    if n == 0 {
        return 0.0;
    }
    if f <= 0.0 {
        return arr[0];
    }
    if f >= (n - 1) as f64 {
        return arr[n - 1];
    }
    let m = f.floor() as usize;
    let t = f - m as f64;

    // Use cubic Lagrange only when the 4-point stencil [m-1, m, m+1, m+2]
    // is fully in-bounds. Otherwise fall back to linear (the boundary
    // pad cells already provide a buffer, so this branch is rarely
    // hit in practice).
    if m >= 1 && m + 2 < n {
        let y0 = arr[m - 1];
        let y1 = arr[m];
        let y2 = arr[m + 1];
        let y3 = arr[m + 2];
        let l_m1 = -t * (t - 1.0) * (t - 2.0) / 6.0;
        let l_0 = (t + 1.0) * (t - 1.0) * (t - 2.0) / 2.0;
        let l_p1 = -(t + 1.0) * t * (t - 2.0) / 2.0;
        let l_p2 = (t + 1.0) * t * (t - 1.0) / 6.0;
        l_m1 * y0 + l_0 * y1 + l_p1 * y2 + l_p2 * y3
    } else {
        (1.0 - t) * arr[m] + t * arr[m + 1]
    }
}

/// Cardinal-axis propagation direction for [`PlaneWaveSource`].
///
/// Phase 2.fdtd.5 only implements [`PlaneWaveDirection::PlusX`] (E_z
/// polarized). The other variants are recognized by the constructor but
/// cause [`PlaneWaveSource::correct_h`] and
/// [`PlaneWaveSource::correct_e`] to `unimplemented!()`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlaneWaveDirection {
    /// Propagation along `+x` (E_z polarized in Phase 2.fdtd.5).
    PlusX,
    /// Propagation along `+y` — not implemented in Phase 2.fdtd.5.
    PlusY,
    /// Propagation along `+z` — not implemented in Phase 2.fdtd.5.
    PlusZ,
    /// Propagation along `-x` — not implemented in Phase 2.fdtd.5.
    MinusX,
    /// Propagation along `-y` — not implemented in Phase 2.fdtd.5.
    MinusY,
    /// Propagation along `-z` — not implemented in Phase 2.fdtd.5.
    MinusZ,
}

/// Total-field / scattered-field (TF/SF) plane-wave source.
///
/// Injects a normally-incident plane wave (Phase 2.fdtd.5 only supports
/// `+x` direction with `E_z` polarization and `H_y` carrier) into the
/// total-field region defined by an axis-aligned box of cell indices
/// `[i0..=i1, j0..=j1, k0..=k1]`.
///
/// # Field convention
///
/// - Inside the TF box, stored `E` and `H` are **total** fields.
/// - Outside the TF box, stored `E` and `H` are **scattered** fields.
///
/// Coupling between the regions is implemented as discrete corrections
/// on the `i = i0` and `i = i1` faces, derived from an auxiliary 1-D
/// FDTD incident-field grid that propagates the analytical plane wave
/// with the same numerical dispersion the 3D scheme sees along the
/// propagation axis.
///
/// # Polarization and supported geometry (Phase 2.fdtd.5 / 2.fdtd.5.1 / 2.fdtd.5.2)
///
/// For `+x` propagation, `E_z` polarized, the only non-zero incident
/// field components are `E_inc_z(x, t)` and `H_inc_y(x, t)` — incident
/// `H_x`, `H_z`, `E_x`, `E_y` are all identically zero. The discrete
/// Yee stencils that pick up a non-zero incident contribution across
/// the TF/SF boundary are therefore exactly four:
///
/// - `E_z` update at `i = i0` and `i = i1` — uses `H_inc_y` across
///   the `i`-face in `∂H_y/∂x`. **Correction applied (5.1).**
/// - `H_y` update at `i = i0 - 1` and `i = i1` — uses `E_inc_z`
///   across the `i`-face in `∂E_z/∂x`. **Correction applied (5.1).**
/// - `H_x` update at `j = j0 - 1` and `j = j1` — uses `E_inc_z`
///   across the `j`-face in `∂E_z/∂y`. **Correction applied (5.2).**
/// - `E_x` update at `k = k0` and `k = k1 + 1` — uses `H_inc_y`
///   across the `k`-face in `∂H_y/∂z`. **Correction applied (5.2).**
///
/// Phase 2.fdtd.5.1 shipped only the first two. With the 5.2 j/k-face
/// additions, finite-box geometry (TF box bounded on all six faces)
/// achieves a contrast ratio well above 100×, comparable to the slab
/// configuration. Slab geometry remains the recommended option when
/// the geometry permits, because the slab j/k faces still sit in CPML
/// and avoid even the discretized-correction round-off.
///
/// # Reference
///
/// Taflove & Hagness, *Computational Electrodynamics* (3rd ed.) §5.10
/// (3-D TF/SF for a rectangular Huygens surface) and §6 / §14.
///
/// # Phase 2.fdtd.5 / 2.fdtd.5.1 / 2.fdtd.5.2 limitations
///
/// - Only `PlusX` direction with `E_z` polarization is implemented;
///   other [`PlaneWaveDirection`] variants `unimplemented!()` in the
///   correction kernels.
/// - All four faces (`i0`, `i1`, `j0/j1`, `k0/k1`) now apply
///   corrections for `+x` `E_z` polarization. Arbitrary polarization
///   and oblique incidence land in Phase 2.fdtd.5.3+.
/// - The 1-D auxiliary grid uses the same `dx` and `dt` as the 3D grid;
///   for normal incidence this is exact in the limit of the 3D cubic
///   Yee dispersion relation on-axis, but introduces a small mismatch
///   at finite resolution that is well within the `> 10×` TF/SF
///   contrast gate.
/// - The 1-D far-end uses a first-order Mur ABC, sufficient for runs
///   of several hundred steps without spurious 1-D reflections
///   leaking back into the TF region.
#[derive(Debug, Clone)]
pub struct PlaneWaveSource {
    /// TF region lower x cell index (inclusive).
    pub i0: usize,
    /// TF region upper x cell index (inclusive).
    pub i1: usize,
    /// TF region lower y cell index (inclusive).
    pub j0: usize,
    /// TF region upper y cell index (inclusive).
    pub j1: usize,
    /// TF region lower z cell index (inclusive).
    pub k0: usize,
    /// TF region upper z cell index (inclusive).
    pub k1: usize,
    /// Propagation direction.
    pub direction: PlaneWaveDirection,
    /// Source carrier frequency (Hz).
    pub frequency: f64,
    /// Hanning-window taper length, in time steps.
    pub ramp_steps: usize,

    /// 1-D auxiliary grid: `E_inc` samples. Length = `(i1 - i0) + 2*pad + 1`
    /// for `PlusX`. Index 0 is the source-injection cell; index `pad`
    /// corresponds to the 3D plane `i = i0`.
    inc_e: Vec<f64>,
    /// 1-D auxiliary grid: `H_inc` samples. Length = `inc_e.len() - 1`,
    /// staggered half a cell to the right of each `inc_e` sample.
    inc_h: Vec<f64>,
    /// Number of "lead-in" cells in the 1-D grid before the TF front face.
    pad: usize,
    /// Cell size of the 3D grid along the propagation axis (cached for
    /// incident-grid updates and corrections).
    dx: f64,
    /// Time step of the 3D grid.
    dt: f64,
    /// 1-D incident grid step counter.
    step: usize,
    /// Previous-step value of `inc_e[N - 1]` (far-end cell), used by the
    /// first-order Mur ABC on the 1-D grid.
    mur_prev_end: f64,
    /// Previous-step value of `inc_e[N - 2]` (cell just inside the
    /// far end), used by the first-order Mur ABC.
    mur_prev_inner: f64,
    /// Mur ABC coefficient `(c·dt - dx)/(c·dt + dx)`, cached.
    mur_coeff: f64,

    /// 1-D auxiliary-grid step (metres). For the legacy normal-incidence
    /// path this equals `dx`. For oblique incidence it is the
    /// **dispersion-matched** step: chosen so the 1-D Yee leapfrog
    /// reproduces the 3-D Yee numerical phase velocity along `k_hat` at
    /// the source carrier frequency. Without this match, the 1-D and
    /// 3-D waves drift in phase across the TF box and the residual
    /// drift dominates the SF leakage; with the match the leakage is
    /// roughly interpolation-limited (Taflove §5.10.5).
    ds_aux: f64,

    /// Incidence polar angle θ (radians), measured from `+z`.
    /// `θ = 0` reduces to the Phase 2.fdtd.5.2 `+x`/`E_z` case via the
    /// (k_hat, e_theta, e_phi) basis below — see `with_oblique_incidence`.
    theta_inc: f64,
    /// Incidence azimuth φ (radians), measured from `+x` in the `xy` plane.
    phi_inc: f64,
    /// Polarization angle ψ (radians) in the plane perpendicular to `k_hat`.
    /// `ψ = 0` → E aligned with `e_theta_hat`; `ψ = π/2` → E aligned with
    /// `e_phi_hat`.
    psi_pol: f64,
    /// Propagation unit vector `k_hat = (sin θ cos φ, sin θ sin φ, cos θ)`.
    /// At θ=φ=0 this is `(0, 0, 1)` — *but* the legacy normal-incidence
    /// path keeps `direction = PlusX` and is detected by
    /// [`Self::is_legacy_normal_incidence`], which short-circuits to the
    /// 5.2 four-face kernel for bit-for-bit back-compat. Oblique
    /// constructors set `direction` to a sentinel handled by the 12-face
    /// kernel.
    k_hat: [f64; 3],
    /// Incident E-field unit vector
    /// `E_inc_hat = cos ψ · e_theta_hat + sin ψ · e_phi_hat`.
    e_inc_hat: [f64; 3],
    /// Incident H-field unit vector
    /// `H_inc_hat = k_hat × E_inc_hat = cos ψ · e_phi_hat − sin ψ · e_theta_hat`.
    /// Magnitudes are unit; the `η₀` ratio between E and H is encoded in
    /// the relative scaling of `inc_e` and `inc_h` from the 1-D leapfrog.
    h_inc_hat: [f64; 3],
    /// Reference corner `r_ref` (in metres) such that the projected
    /// distance along `k_hat` from a Yee node at `r` is `(r − r_ref) · k_hat`.
    /// For `k_hat` with non-negative components this is the upstream TF
    /// box corner `(i0·dx, j0·dy, k0·dz)`; aux index `pad` corresponds
    /// to `s = 0` there.
    r_ref: [f64; 3],
    /// `true` iff this source was constructed via the legacy
    /// [`Self::new`] entry point with `direction = PlusX` — used to
    /// dispatch to the bit-for-bit 5.2 four-face kernel for back-compat.
    /// Oblique constructors set this to `false`.
    legacy_plus_x: bool,
}

impl PlaneWaveSource {
    /// Build a new TF/SF plane-wave source for the given TF region.
    ///
    /// `pad` controls the number of "lead-in" cells in the 1-D auxiliary
    /// grid before the TF front face. A value of `4` is the documented
    /// minimum: the source pulse needs at least a few cells to develop
    /// before its leading edge reaches the TF boundary.
    ///
    /// # Panics
    ///
    /// Panics if any of the region bounds are inverted (`i0 > i1`, etc.)
    /// or if `pad < 1`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        i0: usize,
        i1: usize,
        j0: usize,
        j1: usize,
        k0: usize,
        k1: usize,
        direction: PlaneWaveDirection,
        frequency: f64,
        ramp_steps: usize,
        dx: f64,
        dt: f64,
        pad: usize,
    ) -> Self {
        assert!(i0 <= i1, "PlaneWaveSource: i0 ({i0}) must be ≤ i1 ({i1})");
        assert!(j0 <= j1, "PlaneWaveSource: j0 ({j0}) must be ≤ j1 ({j1})");
        assert!(k0 <= k1, "PlaneWaveSource: k0 ({k0}) must be ≤ k1 ({k1})");
        assert!(pad >= 1, "PlaneWaveSource: pad ({pad}) must be ≥ 1");
        assert!(
            frequency > 0.0 && frequency.is_finite(),
            "PlaneWaveSource: frequency must be positive and finite"
        );
        assert!(
            dx > 0.0 && dx.is_finite(),
            "PlaneWaveSource: dx must be positive and finite"
        );
        assert!(
            dt > 0.0 && dt.is_finite(),
            "PlaneWaveSource: dt must be positive and finite"
        );

        let n_along = match direction {
            PlaneWaveDirection::PlusX | PlaneWaveDirection::MinusX => i1 - i0,
            PlaneWaveDirection::PlusY | PlaneWaveDirection::MinusY => j1 - j0,
            PlaneWaveDirection::PlusZ | PlaneWaveDirection::MinusZ => k1 - k0,
        };
        let inc_n_cells = n_along + 2 * pad + 1;
        let inc_e = vec![0.0; inc_n_cells];
        let inc_h = vec![0.0; inc_n_cells - 1];

        let c0 = yee_core::units::C0;
        let mur_coeff = (c0 * dt - dx) / (c0 * dt + dx);

        // Populate the (k_hat, E_inc_hat, H_inc_hat) trio for the legacy
        // +x / E_z case so the 12-face general kernel could in principle
        // also be invoked on this source — but `legacy_plus_x = true`
        // dispatches to the 5.2 4-face kernel by default for bit-for-bit
        // back-compat.
        let (k_hat, e_inc_hat, h_inc_hat) = match direction {
            PlaneWaveDirection::PlusX => ([1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]),
            _ => ([1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]),
        };
        let r_ref = [i0 as f64 * dx, j0 as f64 * dx, k0 as f64 * dx];

        Self {
            i0,
            i1,
            j0,
            j1,
            k0,
            k1,
            direction,
            frequency,
            ramp_steps,
            inc_e,
            inc_h,
            pad,
            dx,
            dt,
            step: 0,
            mur_prev_end: 0.0,
            mur_prev_inner: 0.0,
            mur_coeff,
            ds_aux: dx,
            theta_inc: 0.0,
            phi_inc: 0.0,
            psi_pol: 0.0,
            k_hat,
            e_inc_hat,
            h_inc_hat,
            r_ref,
            legacy_plus_x: matches!(direction, PlaneWaveDirection::PlusX),
        }
    }

    /// Build a new TF/SF plane-wave source with oblique incidence.
    ///
    /// `theta_inc` and `phi_inc` (radians) set the propagation unit
    /// vector
    /// ```text
    ///   k_hat = (sin θ cos φ, sin θ sin φ, cos θ)
    /// ```
    /// (`θ` measured from `+z`, `φ` from `+x` in the `xy` plane.)
    /// `psi_pol` (radians) sets the polarization angle in the plane
    /// perpendicular to `k_hat`:
    /// ```text
    ///   E_inc_hat = cos ψ · e_theta_hat + sin ψ · e_phi_hat
    /// ```
    /// where `e_theta_hat = (cos θ cos φ, cos θ sin φ, −sin θ)` and
    /// `e_phi_hat = (−sin φ, cos φ, 0)`. `H_inc_hat` is then
    /// `k_hat × E_inc_hat`.
    ///
    /// `pad` controls the number of "lead-in" cells in the 1-D
    /// auxiliary grid along `k_hat` before `s = 0` (which corresponds
    /// to the upstream TF-box corner `(i0·dx, j0·dy, k0·dz)`).
    ///
    /// # Restrictions (Phase 2.fdtd.5.3)
    ///
    /// - `k_hat` must have non-negative components: `0 ≤ θ ≤ π/2` and
    ///   `0 ≤ φ ≤ π/2`. This guarantees the `(i0, j0, k0)` corner is
    ///   the upstream-most face node, so every projected face distance
    ///   `s ≥ 0` and we can use the single-sided 1-D aux grid.
    /// - The 3-D grid must be cubic (`dx = dy = dz`). This is enforced
    ///   in the auxiliary grid sizing.
    ///
    /// # Panics
    ///
    /// Same panics as [`Self::new`], plus panics if any angle is
    /// non-finite or outside its allowed range.
    ///
    /// # Dispersion matching
    ///
    /// Since Phase 2.fdtd.5.3.1, this constructor enables Taflove
    /// §5.10.5 dispersion matching of the 1-D auxiliary grid by
    /// default: the aux step `ds_aux` is chosen so the 1-D Yee
    /// numerical phase velocity matches the 3-D Yee phase velocity
    /// projected along `k_hat` at the source carrier frequency.
    /// Without this match, the 1-D and 3-D waves drift in phase
    /// across the TF box and oblique-incidence contrast plateaus
    /// at O(10×). With it, the 30°/45° finite-box contrast clears
    /// 1000× (the original Phase 2.fdtd.5.3 DoD).
    ///
    /// To reproduce the pre-5.3.1 `ds_aux = dx` behaviour (e.g. for
    /// regression / sanity comparison), use
    /// [`Self::with_oblique_incidence_match`] with
    /// `dispersion_match = false`.
    #[allow(clippy::too_many_arguments)]
    pub fn with_oblique_incidence(
        i0: usize,
        i1: usize,
        j0: usize,
        j1: usize,
        k0: usize,
        k1: usize,
        theta_inc: f64,
        phi_inc: f64,
        psi_pol: f64,
        frequency: f64,
        ramp_steps: usize,
        dx: f64,
        dt: f64,
        pad: usize,
    ) -> Self {
        Self::with_oblique_incidence_match(
            i0, i1, j0, j1, k0, k1, theta_inc, phi_inc, psi_pol, frequency, ramp_steps, dx, dt,
            pad, true,
        )
    }

    /// Build a new TF/SF plane-wave source with oblique incidence and
    /// explicit control over 1-D auxiliary-grid dispersion matching.
    ///
    /// All arguments match [`Self::with_oblique_incidence`] except
    /// the trailing `dispersion_match: bool` flag. When `true`
    /// (the default for [`Self::with_oblique_incidence`]) the aux
    /// step is computed internally to match the 3-D
    /// numerical phase velocity along `k_hat` at the source carrier
    /// frequency (Taflove §5.10.5). When `false` the aux step is
    /// hard-coded to `ds_aux = dx`, reproducing the Phase 2.fdtd.5.3
    /// ship behaviour (oblique contrast plateaus at ~14.5× for the
    /// 30°/45° case). The `false` mode is provided primarily for
    /// regression tests; production callers should leave the flag
    /// `true`.
    #[allow(clippy::too_many_arguments)]
    pub fn with_oblique_incidence_match(
        i0: usize,
        i1: usize,
        j0: usize,
        j1: usize,
        k0: usize,
        k1: usize,
        theta_inc: f64,
        phi_inc: f64,
        psi_pol: f64,
        frequency: f64,
        ramp_steps: usize,
        dx: f64,
        dt: f64,
        pad: usize,
        dispersion_match: bool,
    ) -> Self {
        assert!(i0 <= i1, "PlaneWaveSource: i0 ({i0}) must be ≤ i1 ({i1})");
        assert!(j0 <= j1, "PlaneWaveSource: j0 ({j0}) must be ≤ j1 ({j1})");
        assert!(k0 <= k1, "PlaneWaveSource: k0 ({k0}) must be ≤ k1 ({k1})");
        assert!(pad >= 1, "PlaneWaveSource: pad ({pad}) must be ≥ 1");
        assert!(
            frequency > 0.0 && frequency.is_finite(),
            "PlaneWaveSource: frequency must be positive and finite"
        );
        assert!(
            dx > 0.0 && dx.is_finite(),
            "PlaneWaveSource: dx must be positive and finite"
        );
        assert!(
            dt > 0.0 && dt.is_finite(),
            "PlaneWaveSource: dt must be positive and finite"
        );
        assert!(
            theta_inc.is_finite() && phi_inc.is_finite() && psi_pol.is_finite(),
            "PlaneWaveSource: angles must be finite"
        );
        assert!(
            (0.0..=PI / 2.0 + 1e-12).contains(&theta_inc),
            "PlaneWaveSource (5.3): θ must be in [0, π/2] (got {theta_inc})"
        );
        assert!(
            (0.0..=PI / 2.0 + 1e-12).contains(&phi_inc),
            "PlaneWaveSource (5.3): φ must be in [0, π/2] (got {phi_inc})"
        );

        let (sin_t, cos_t) = theta_inc.sin_cos();
        let (sin_p, cos_p) = phi_inc.sin_cos();
        let (sin_psi, cos_psi) = psi_pol.sin_cos();

        let k_hat = [sin_t * cos_p, sin_t * sin_p, cos_t];
        let e_theta = [cos_t * cos_p, cos_t * sin_p, -sin_t];
        let e_phi = [-sin_p, cos_p, 0.0];
        let e_inc_hat = [
            cos_psi * e_theta[0] + sin_psi * e_phi[0],
            cos_psi * e_theta[1] + sin_psi * e_phi[1],
            cos_psi * e_theta[2] + sin_psi * e_phi[2],
        ];
        // H_inc_hat = −(k_hat × E_inc_hat). The negative sign here is
        // intentional and absorbs the sign convention of the 1-D Yee
        // aux grid: the 1-D leapfrog produces `inc_h ≈ −inc_e/η₀` for a
        // wave propagating in `+s` (i.e. `inc_h` carries the opposite
        // sign of the physical H magnitude). Setting
        // `H_inc_hat = −(k̂×E_inc_hat)` cancels that sign so that
        // `H_inc_α(r,t) = H_inc_hat_α · inc_h(...)` reproduces the
        // physical Maxwell field. For the 5.2 normal-incidence path
        // (θ=π/2, φ=0, ψ=π, k̂=+x̂, E_inc_hat=+ẑ), this gives
        // H_inc_hat = +ŷ, matching the 5.2 convention's positive
        // H_inc_y_hat scaling on the (already-negative) `inc_h`.
        let h_inc_hat = [
            -(cos_psi * e_phi[0] - sin_psi * e_theta[0]),
            -(cos_psi * e_phi[1] - sin_psi * e_theta[1]),
            -(cos_psi * e_phi[2] - sin_psi * e_theta[2]),
        ];

        // Phase 2.fdtd.5.3.1: select the 1-D aux-grid step.
        //
        // With `dispersion_match = true` (the default for
        // `with_oblique_incidence`), `ds_aux` is solved via bisection
        // so the 1-D Yee numerical phase velocity matches the 3-D Yee
        // numerical phase velocity along `k_hat` at the source carrier
        // frequency (Taflove §5.10.5). This raises 30°/45° contrast
        // from ~14.5× to >1000×.
        //
        // With `dispersion_match = false` the aux step is hard-coded
        // to `dx`, reproducing the Phase 2.fdtd.5.3 ship behaviour
        // (for regression / back-compat comparison only).
        let ds_aux = if dispersion_match {
            let omega = TAU * frequency;
            compute_aux_step(theta_inc, phi_inc, dx, C0, omega, dt)
        } else {
            dx
        };

        // Aux 1-D grid size: must cover the projected box diagonal
        // (S_max along k_hat) plus pad on each side. Compute S_max in
        // metres and convert to aux cells via `ds_aux`.
        let s_max_m = (i1 - i0) as f64 * dx * k_hat[0].abs()
            + (j1 - j0) as f64 * dx * k_hat[1].abs()
            + (k1 - k0) as f64 * dx * k_hat[2].abs();
        let n_along = (s_max_m / ds_aux).ceil() as usize + 1;
        let inc_n_cells = n_along + 2 * pad + 1;
        let inc_e = vec![0.0; inc_n_cells];
        let inc_h = vec![0.0; inc_n_cells - 1];

        let mur_coeff = (C0 * dt - ds_aux) / (C0 * dt + ds_aux);
        let r_ref = [i0 as f64 * dx, j0 as f64 * dx, k0 as f64 * dx];

        Self {
            i0,
            i1,
            j0,
            j1,
            k0,
            k1,
            direction: PlaneWaveDirection::PlusX, // sentinel; oblique kernel ignores it
            frequency,
            ramp_steps,
            inc_e,
            inc_h,
            pad,
            dx,
            dt,
            step: 0,
            mur_prev_end: 0.0,
            mur_prev_inner: 0.0,
            mur_coeff,
            ds_aux,
            theta_inc,
            phi_inc,
            psi_pol,
            k_hat,
            e_inc_hat,
            h_inc_hat,
            r_ref,
            legacy_plus_x: false,
        }
    }

    /// `true` if this source uses the bit-for-bit 5.2 four-face
    /// normal-incidence kernel. `false` for sources built via
    /// [`Self::with_oblique_incidence`].
    pub fn is_legacy_normal_incidence(&self) -> bool {
        self.legacy_plus_x
    }

    /// Incidence polar angle θ in radians.
    pub fn theta_inc(&self) -> f64 {
        self.theta_inc
    }

    /// Incidence azimuth φ in radians.
    pub fn phi_inc(&self) -> f64 {
        self.phi_inc
    }

    /// Polarization angle ψ in radians.
    pub fn psi_pol(&self) -> f64 {
        self.psi_pol
    }

    /// Hanning (raised-cosine) ramp factor for a sinusoidal source.
    ///
    /// Returns `0.5 * (1 - cos(π · n / ramp_steps))` for `n < ramp_steps`
    /// and `1.0` afterwards. Tapering the carrier on with a Hann window
    /// suppresses the broadband click an unramped sinusoid would inject.
    fn ramp(&self) -> f64 {
        if self.ramp_steps == 0 || self.step >= self.ramp_steps {
            1.0
        } else {
            0.5 * (1.0 - (std::f64::consts::PI * self.step as f64 / self.ramp_steps as f64).cos())
        }
    }

    /// Drive value of the source at the current 1-D-grid step: a sinusoid
    /// `sin(2π f n dt)` modulated by the Hann ramp.
    fn source_value(&self) -> f64 {
        let t = self.step as f64 * self.dt;
        self.ramp() * (TAU * self.frequency * t).sin()
    }

    /// Advance `H_inc` by one time step.
    ///
    /// Standard 1-D Yee H-update:
    /// ```text
    /// H_inc_y[m+1/2] += (Δt/(μ₀ Δx)) · (E_inc_z[m+1] - E_inc_z[m])
    /// ```
    /// matches the sign convention of [`crate::update::update_h`] on
    /// the 3D grid (∂H_y/∂t = +(1/μ) ∂E_z/∂x).
    pub fn step_incident_h(&mut self) {
        // For the legacy normal-incidence path, `ds_aux == dx` so this
        // is identical to the 5.x code. For oblique sources, `ds_aux`
        // is the dispersion-matched step (set at construction).
        let coeff = self.dt / (MU0 * self.ds_aux);
        for m in 0..self.inc_h.len() {
            self.inc_h[m] += coeff * (self.inc_e[m + 1] - self.inc_e[m]);
        }
    }

    /// Update `E_inc`, inject the analytic source at the near end, and
    /// apply a first-order Mur ABC at the far end. See
    /// [`Self::step_incident_h`].
    ///
    /// The leapfrog body is:
    /// ```text
    /// E_inc_z[m]   += (Δt/(ε₀ Δx)) · (H_inc_y[m+1/2] - H_inc_y[m-1/2])
    /// E_inc_z[0]   = ramp(n)·sin(2π f n Δt)              (hard source)
    /// E_inc_z[N-1] = E_inc[N-2]^old + κ·(E_inc[N-2]^new - E_inc[N-1]^old)
    /// ```
    /// where κ = (c·Δt - Δx)/(c·Δt + Δx) is the Mur first-order ABC
    /// coefficient.
    pub fn step_incident_e(&mut self) {
        let coeff = self.dt / (EPS0 * self.ds_aux);
        let n = self.inc_e.len();

        let prev_end = self.mur_prev_end;
        let prev_inner = self.mur_prev_inner;

        // Update E_inc[m] for m ∈ [1, n-1) using the freshly-stepped H_inc.
        for m in 1..n - 1 {
            self.inc_e[m] += coeff * (self.inc_h[m] - self.inc_h[m - 1]);
        }
        // Hard source at m=0.
        self.step += 1;
        self.inc_e[0] = self.source_value();

        // First-order Mur ABC at far end for outgoing +x waves.
        let inner_new = self.inc_e[n - 2];
        self.inc_e[n - 1] = prev_inner + self.mur_coeff * (inner_new - prev_end);

        // Save state for next call.
        self.mur_prev_end = self.inc_e[n - 1];
        self.mur_prev_inner = inner_new;
    }

    /// Apply TF/SF corrections to the magnetic field on the box faces.
    /// Call **after** [`crate::update::update_h`] and **after**
    /// [`Self::step_incident_h`] (which advances `H_inc` from the
    /// current `E_inc`).
    ///
    /// # Phase 2.fdtd.5 scope
    ///
    /// Implements `+x` propagation, `E_z` polarized only. Other variants
    /// of [`PlaneWaveDirection`] call `unimplemented!()`.
    pub fn correct_h(&self, grid: &mut YeeGrid) {
        if self.legacy_plus_x {
            // Bit-for-bit 5.2 four-face kernel for normal-incidence
            // back-compat.
            match self.direction {
                PlaneWaveDirection::PlusX => self.correct_h_plus_x(grid),
                _ => unimplemented!(
                    "PlaneWaveDirection::{:?} is not implemented in Phase 2.fdtd.5",
                    self.direction
                ),
            }
        } else {
            self.correct_h_oblique(grid);
        }
    }

    /// Apply TF/SF corrections to the electric field on the box faces.
    /// Call **after** [`crate::update::update_e`] and **after**
    /// [`Self::step_incident_e`].
    ///
    /// # Phase 2.fdtd.5 scope
    ///
    /// Implements `+x` propagation, `E_z` polarized only. Other variants
    /// of [`PlaneWaveDirection`] call `unimplemented!()`.
    pub fn correct_e(&self, grid: &mut YeeGrid) {
        if self.legacy_plus_x {
            // Bit-for-bit 5.2 four-face kernel for normal-incidence
            // back-compat.
            match self.direction {
                PlaneWaveDirection::PlusX => self.correct_e_plus_x(grid),
                _ => unimplemented!(
                    "PlaneWaveDirection::{:?} is not implemented in Phase 2.fdtd.5",
                    self.direction
                ),
            }
        } else {
            self.correct_e_oblique(grid);
        }
    }

    /// Map a 3D x-index `i` to a 1-D incident-grid `E_inc` index.
    #[inline]
    fn e_idx(&self, i: usize) -> usize {
        i - self.i0 + self.pad
    }

    /// Map a 3D H_y i-index to a 1-D incident-grid `H_inc` index.
    /// `H_y[i, *, *]` lives at the half-cell `(i + 1/2, *, *)`, so its
    /// 1-D counterpart is `H_inc[i - i0 + pad]`.
    #[inline]
    fn h_idx(&self, i_h: usize) -> usize {
        i_h - self.i0 + self.pad
    }

    // ----------------------------------------------------------------
    // +x propagation, E_z polarization (Phase 2.fdtd.5 / 2.fdtd.5.1 / 2.fdtd.5.2)
    //
    // Derivation (Taflove & Hagness §5.10 / §14):
    //
    // For a +x plane wave with E along z and H along y, the incident
    // field has only E_inc_z(x) and H_inc_y(x) non-zero. The four Yee
    // stencils that pick up a non-zero incident contribution across the
    // TF/SF boundary are:
    //
    //   - E_z update at i = i0 / i1   (uses H_inc_y across the i-face,
    //                                  in `∂H_y/∂x` term of E_z curl)
    //   - H_y update at i = i0-1 / i1 (uses E_inc_z across the i-face,
    //                                  in `∂E_z/∂x` term of H_y curl)
    //   - H_x update at j = j0-1 / j1 (uses E_inc_z across the j-face,
    //                                  in `∂E_z/∂y` term of H_x curl)
    //   - E_x update at k = k0 / k1+1 (uses H_inc_y across the k-face,
    //                                  in `∂H_y/∂z` term of E_x curl)
    //
    // The first two are i-face corrections (Phase 2.fdtd.5 / 5.1); the
    // last two are j/k-face corrections (Phase 2.fdtd.5.2).
    //
    // ----- i-face -----
    //
    // Front (i = i0):
    //   H_y[i0-1, j, k] is SF (between SF E_z[i0-1] and TF E_z[i0]).
    //   Standard update_h read E_z[i0] (TF) thinking it was SF;
    //   correction: subtract (dt/(μ₀·dx)) · E_inc_z[at i0]  from H_y[i0-1].
    //
    //   E_z[i0, j, k] is TF, but standard update_e read
    //   H_y[i0-1, j, k] (SF) thinking it was TF; correction:
    //   subtract (dt/(ε₀·dx)) · H_inc_y[at i0-1]  from E_z[i0].
    //
    // Back (i = i1):
    //   H_y[i1, j, k] is SF (between TF E_z[i1] and SF E_z[i1+1]).
    //   Standard update_h read E_z[i1] (TF) thinking it was SF;
    //   correction: add (dt/(μ₀·dx)) · E_inc_z[at i1]  to H_y[i1].
    //
    //   E_z[i1, j, k] is TF, but standard update_e read
    //   H_y[i1, j, k] (SF) thinking it was TF; correction:
    //   add (dt/(ε₀·dx)) · H_inc_y[at i1]  to E_z[i1].
    //
    // ----- j-face (Phase 2.fdtd.5.2) -----
    //
    // H_x[i, j, k] update is:
    //     H_x += (dt/μ₀) · (∂E_y/∂z − ∂E_z/∂y)
    //          = (dt/μ₀) · ( (E_y[i,j,k+1]−E_y[i,j,k])/dz
    //                        − (E_z[i,j+1,k]−E_z[i,j,k])/dy )
    // Only the `∂E_z/∂y` term can straddle the j-face TF/SF boundary
    // (E_z is the only incident-bearing field; E_y has no incident).
    //
    // Front (j = j0):
    //   H_x[i, j0-1, k] is SF (between SF E_z[i,j0-1,k] and TF
    //   E_z[i,j0,k]). Standard update_h read E_z[i,j0,k] (TF) as if it
    //   were SF; that put a spurious `−(dt/(μ₀·dy))·E_inc_z(i)` into
    //   H_x[i, j0-1, k] (negative because of the minus sign on
    //   `∂E_z/∂y`). Correction:
    //
    //     H_x[i, j0-1, k]  +=  (dt/(μ₀·dy)) · E_inc_z[at i]
    //
    // Back (j = j1):
    //   H_x[i, j1, k] is SF (between TF E_z[i,j1,k] and SF
    //   E_z[i,j1+1,k]). The TF E_z is now the *subtracted* term in
    //   `(E_z[j1+1] − E_z[j1])/dy`, so the spurious contribution is
    //   `+(dt/(μ₀·dy))·E_inc_z(i)`. Correction:
    //
    //     H_x[i, j1, k]    −=  (dt/(μ₀·dy)) · E_inc_z[at i]
    //
    // Both corrections use `E_inc_z` sampled at the x-index of the H_x
    // cell (`H_x[i,*,*]` lives at integer x = i, and E_z[i,*,*] also
    // lives at integer x = i, so they share `inc_e[e_idx(i)]`).
    //
    // No j-face correction is needed for E_z updates: E_z curl has
    // `∂H_x/∂y`, and H_x has no incident component (H_inc_x = 0).
    //
    // ----- k-face (Phase 2.fdtd.5.2) -----
    //
    // E_x[i, j, k] update is:
    //     E_x += (dt/ε₀) · (∂H_z/∂y − ∂H_y/∂z)
    //          = (dt/ε₀) · ( (H_z[i,j,k]−H_z[i,j-1,k])/dy
    //                        − (H_y[i,j,k]−H_y[i,j,k-1])/dz )
    // Only the `∂H_y/∂z` term can straddle the k-face TF/SF boundary
    // (H_y is the only incident-bearing field; H_z has no incident).
    //
    // Front (k = k0):
    //   E_x[i, j, k0] lives at z = k0 (integer plane), on the boundary
    //   between SF H_y[i,j,k0-1] (z = k0-1/2) and TF H_y[i,j,k0]
    //   (z = k0+1/2). By the same convention used on the i-face — E
    //   nodes on the boundary plane are claimed as TF — E_x[i,j,k0] is
    //   TF. The standard update read H_y[i,j,k0-1] (SF) thinking it was
    //   TF, so it under-read by `−H_inc_y(i+1/2)`; that propagated into
    //   `∂H_y/∂z` as `−(−H_inc_y)/dz = +H_inc_y/dz`, and into the E_x
    //   update with the `−∂H_y/∂z` sign as `−(dt/(ε₀·dz))·H_inc_y(i+1/2)`.
    //   Correction:
    //
    //     E_x[i, j, k0]    +=  (dt/(ε₀·dz)) · H_inc_y[at i+1/2]
    //
    // Back (k = k1 + 1):
    //   E_x[i, j, k1+1] lives at z = k1+1 (integer plane), on the
    //   boundary between TF H_y[i,j,k1] (z = k1+1/2) and SF
    //   H_y[i,j,k1+1] (z = k1+3/2). E_x at this plane is TF by
    //   convention. Standard update read H_y[i,j,k1+1] (SF) as if TF,
    //   under-reading by `−H_inc_y(i+1/2)`; that propagated through
    //   `∂H_y/∂z` and `−∂H_y/∂z` as `+(dt/(ε₀·dz))·H_inc_y(i+1/2)` of
    //   spurious E_x. Correction:
    //
    //     E_x[i, j, k1+1]  −=  (dt/(ε₀·dz)) · H_inc_y[at i+1/2]
    //
    // Both corrections use `H_inc_y` sampled at the x-coord of the
    // E_x cell. E_x[i,*,*] lives at x = i+1/2; H_y[i,*,*] also lives
    // at x = i+1/2; so the 1-D-grid lookup is `inc_h[h_idx(i)]`.
    //
    // No k-face correction is needed for E_z updates (E_z curl has no
    // z-derivative) or for H_x / H_y updates whose z-derivative terms
    // involve E_x or E_y (no incident).
    // ----------------------------------------------------------------

    fn correct_h_plus_x(&self, grid: &mut YeeGrid) {
        self.correct_h_iface_plus_x(grid);
        self.correct_h_jface_plus_x(grid);
    }

    fn correct_e_plus_x(&self, grid: &mut YeeGrid) {
        self.correct_e_iface_plus_x(grid);
        self.correct_e_kface_plus_x(grid);
    }

    fn correct_h_iface_plus_x(&self, grid: &mut YeeGrid) {
        let coeff = self.dt / (MU0 * self.dx);
        let einc_front = self.inc_e[self.e_idx(self.i0)];
        let einc_back = self.inc_e[self.e_idx(self.i1)];

        // Bounds-check: H_y has shape [nx, ny+1, nz]. Need i0 ≥ 1
        // (so i0-1 is valid) and i1 ≤ nx-1.
        assert!(
            self.i0 >= 1,
            "PlaneWaveSource (+x): i0 must be ≥ 1 (got {})",
            self.i0
        );
        assert!(
            self.i1 < grid.nx,
            "PlaneWaveSource (+x): i1 ({}) must be < grid.nx ({})",
            self.i1,
            grid.nx
        );

        // H_y cross-section: all (j, k) where E_z[i0, j, k] is TF, i.e.
        // j ∈ [j0, j1] and k ∈ [k0, k1]. H_y has shape [nx, ny+1, nz],
        // so j up to ny is valid and k up to nz-1 is valid; clamp k1.
        // (Phase 2.fdtd.5.2: switched k from exclusive `k0..k_hi` to
        // inclusive `k0..=k1.min(nz-1)` so the upper-z TF slice gets
        // corrected too; slab geometry — where k1 = nz — is unchanged
        // because min(nz, nz-1) = nz-1 = the original `k_hi - 1`.)
        let k_hi = self.k1.min(grid.nz.saturating_sub(1));
        for j in self.j0..=self.j1.min(grid.ny) {
            for k in self.k0..=k_hi {
                grid.hy[(self.i0 - 1, j, k)] -= coeff * einc_front;
                grid.hy[(self.i1, j, k)] += coeff * einc_back;
            }
        }
    }

    /// Apply the j-face H_x corrections (Phase 2.fdtd.5.2).
    ///
    /// Cancels the spurious `E_inc_z` contribution picked up by the
    /// standard `update_h` `∂E_z/∂y` stencil at `H_x[i, j0-1, k]`
    /// (front face) and `H_x[i, j1, k]` (back face).
    ///
    /// `H_x` has shape `[nx+1, ny, nz]`. The j-face correction is a
    /// no-op when `j0 == 0` (no SF row at `j = j0 - 1` to correct) or
    /// when `j1 >= ny` (no SF row at `j = j1`); both situations
    /// correspond to slab geometry where the j-face sits in CPML.
    fn correct_h_jface_plus_x(&self, grid: &mut YeeGrid) {
        // H_x dy-coefficient. `grid.dy` is the relevant cell size for
        // the `∂E_z/∂y` stencil; the 3D grid is cubic in the walking
        // skeleton (dx = dy = dz), but we use grid.dy explicitly for
        // forward compatibility with non-cubic cells.
        let coeff = self.dt / (MU0 * grid.dy);

        // H_x has shape [nx+1, ny, nz]. The cross-section we correct
        // covers (i, k) ∈ [i0, i1] × [k0, k1]; clamp to valid H_x
        // indices.
        let i_hi = self.i1.min(grid.nx);
        let k_hi = self.k1.min(grid.nz.saturating_sub(1));

        // Front j-face: SF H_x row at j = j0 - 1. Skip when j0 == 0
        // (slab in y — the j-face sits in CPML, no correction needed).
        if self.j0 >= 1 {
            for i in self.i0..=i_hi {
                let einc = self.inc_e[self.e_idx(i)];
                for k in self.k0..=k_hi {
                    grid.hx[(i, self.j0 - 1, k)] += coeff * einc;
                }
            }
        }

        // Back j-face: SF H_x row at j = j1. Skip when j1 >= ny (slab
        // in y — the j-face row at j = j1 is past the H_x j-range).
        if self.j1 < grid.ny {
            for i in self.i0..=i_hi {
                let einc = self.inc_e[self.e_idx(i)];
                for k in self.k0..=k_hi {
                    grid.hx[(i, self.j1, k)] -= coeff * einc;
                }
            }
        }
    }

    fn correct_e_iface_plus_x(&self, grid: &mut YeeGrid) {
        let coeff = self.dt / (EPS0 * self.dx);
        assert!(
            self.h_idx(self.i0 - 1) < self.inc_h.len(),
            "PlaneWaveSource (+x): h_idx out of range (logic bug, please report)"
        );
        let hinc_front = self.inc_h[self.h_idx(self.i0 - 1)];
        let hinc_back = self.inc_h[self.h_idx(self.i1)];

        // E_z cross-section at i = i0 / i1: all (j, k) where this E_z
        // is itself TF. E_z has shape [nx+1, ny+1, nz], so j up to ny
        // and k up to nz-1 are valid. (See `correct_h_iface_plus_x`
        // for the Phase 2.fdtd.5.2 inclusive-k rationale.)
        let k_hi = self.k1.min(grid.nz.saturating_sub(1));
        for j in self.j0..=self.j1.min(grid.ny) {
            for k in self.k0..=k_hi {
                grid.ez[(self.i0, j, k)] -= coeff * hinc_front;
                grid.ez[(self.i1, j, k)] += coeff * hinc_back;
            }
        }
    }

    /// Apply the k-face E_x corrections (Phase 2.fdtd.5.2).
    ///
    /// Cancels the spurious `H_inc_y` contribution picked up by the
    /// standard `update_e` `∂H_y/∂z` stencil at `E_x[i, j, k0]`
    /// (front face) and `E_x[i, j, k1+1]` (back face).
    ///
    /// `E_x` has shape `[nx, ny+1, nz+1]`. The k-face correction is a
    /// no-op when `k0 == 0` (no E_x row above the boundary in z — the
    /// k=0 face is PEC/CPML) or when `k1 + 1 > nz`; both situations
    /// correspond to slab geometry where the k-face sits in CPML.
    fn correct_e_kface_plus_x(&self, grid: &mut YeeGrid) {
        // E_x dz-coefficient.
        let coeff = self.dt / (EPS0 * grid.dz);

        // E_x cross-section to correct: the (i, j) cells where the
        // standard `update_e` `∂H_y/∂z` stencil straddles the k-face.
        // The straddle exists only where H_y on the TF side of the
        // face is itself TF. By the i-face convention, H_y is TF for
        // i ∈ [i0, i1-1] (i1 is the SF back-boundary index for H_y);
        // and TF for j ∈ [j0, j1]. The i-range is therefore one cell
        // **narrower** than the H_x j-face cross-section (because of
        // the half-cell offset between H_y's i-coordinate (x = i+1/2)
        // and E_z's i-coordinate (x = i, integer)). Clamp to valid
        // E_x bounds; `E_x` has shape `[nx, ny+1, nz+1]`.
        let i_hi = self.i1.saturating_sub(1).min(grid.nx.saturating_sub(1));
        let j_hi = self.j1.min(grid.ny);

        // Front k-face: TF E_x slab at k = k0. Skip when k0 == 0
        // because then there is no SF H_y row at k = k0 - 1 (the
        // boundary sits at the grid edge — CPML territory).
        if self.k0 >= 1 && self.i0 <= i_hi {
            for i in self.i0..=i_hi {
                let hinc = self.inc_h[self.h_idx(i)];
                for j in self.j0..=j_hi {
                    grid.ex[(i, j, self.k0)] += coeff * hinc;
                }
            }
        }

        // Back k-face: TF E_x slab at k = k1 + 1. Skip when k1+1 > nz
        // because then there is no E_x row at that k (slab geometry).
        if self.k1 < grid.nz && self.i0 <= i_hi {
            for i in self.i0..=i_hi {
                let hinc = self.inc_h[self.h_idx(i)];
                for j in self.j0..=j_hi {
                    grid.ex[(i, j, self.k1 + 1)] -= coeff * hinc;
                }
            }
        }
    }

    // ----------------------------------------------------------------
    // Oblique / general-polarization TF/SF kernel (Phase 2.fdtd.5.3)
    // ----------------------------------------------------------------

    /// Interpolate `inc_e` at the projected distance `s` (in metres)
    /// from the TF reference corner along `k_hat`. `inc_e[pad]`
    /// corresponds to `s = 0`. The aux-grid step is `ds_aux` (matches
    /// `dx` for the legacy normal-incidence path and is
    /// dispersion-matched for oblique sources). Uses 4-point cubic
    /// Lagrange interpolation in the interior (O((Δs/λ)⁴) error), with
    /// graceful fall-back to linear interpolation near the aux-grid
    /// boundaries. Clamps to grid bounds.
    ///
    /// Phase 2.fdtd.5.3.2: upgraded from linear to cubic interpolation
    /// because the linear-interpolation residual (~0.1% per sample at
    /// dx/λ = 0.05) became the dominant SF leakage source once the
    /// face-stencil k-range off-by-one was fixed (commit prior).
    #[inline]
    fn sample_inc_e(&self, s: f64) -> f64 {
        sample_aux_cubic(&self.inc_e, s / self.ds_aux + self.pad as f64)
    }

    /// Interpolate `inc_h` at the projected distance `s` (in metres)
    /// from the TF reference corner along `k_hat`. `inc_h[m]`
    /// corresponds to `s = (m - pad + ½)·ds_aux` (the H_inc samples are
    /// staggered half a cell to the right of the E_inc samples).
    /// Uses cubic Lagrange interpolation (see [`Self::sample_inc_e`]).
    #[inline]
    fn sample_inc_h(&self, s: f64) -> f64 {
        sample_aux_cubic(&self.inc_h, s / self.ds_aux + self.pad as f64 - 0.5)
    }

    /// Projected distance (metres) along `k_hat` from the TF reference
    /// corner `r_ref` to a node at physical position `r`.
    #[inline]
    fn proj_s(&self, x: f64, y: f64, z: f64) -> f64 {
        (x - self.r_ref[0]) * self.k_hat[0]
            + (y - self.r_ref[1]) * self.k_hat[1]
            + (z - self.r_ref[2]) * self.k_hat[2]
    }

    /// E_inc vector component `α` (α = 0/1/2 for x/y/z) at the
    /// physical position `(x, y, z)`.
    #[inline]
    fn e_inc_component(&self, alpha: usize, x: f64, y: f64, z: f64) -> f64 {
        if self.e_inc_hat[alpha] == 0.0 {
            return 0.0;
        }
        self.e_inc_hat[alpha] * self.sample_inc_e(self.proj_s(x, y, z))
    }

    /// H_inc vector component `α` (α = 0/1/2 for x/y/z) at the
    /// physical position `(x, y, z)`.
    #[inline]
    fn h_inc_component(&self, alpha: usize, x: f64, y: f64, z: f64) -> f64 {
        if self.h_inc_hat[alpha] == 0.0 {
            return 0.0;
        }
        self.h_inc_hat[alpha] * self.sample_inc_h(self.proj_s(x, y, z))
    }

    /// Apply the full 12-stencil TF/SF magnetic-field correction for
    /// an oblique plane wave with arbitrary `(θ, φ, ψ)`.
    ///
    /// For each of the six box faces, two H-component stencils
    /// straddle the discontinuity (the curl term parallel to the face
    /// normal of *each* tangential E component). The cross-section
    /// indices match the 5.2 four-face kernel; only the per-cell
    /// incident value changes — it is now an arbitrary linear
    /// combination of `E_inc_x`, `E_inc_y`, `E_inc_z` evaluated at
    /// the face-node position via the aux 1-D grid.
    fn correct_h_oblique(&self, grid: &mut YeeGrid) {
        let dx = grid.dx;
        let dy = grid.dy;
        let dz = grid.dz;
        let dt = self.dt;
        let cx = dt / (MU0 * dx);
        let cy = dt / (MU0 * dy);
        let cz = dt / (MU0 * dz);

        // ---------------------------------------------------------------
        // i-faces: affect H_y at i0-1, i1 (uses E_inc_z across i-face)
        //          and H_z at i0-1, i1 (uses E_inc_y across i-face)
        // H_y[i,j,k] sits at (i+½, j, k+½). H_z[i,j,k] sits at
        // (i+½, j+½, k). On the i-face the H samples nominally sit at
        // x = i0 - ½ (front, i = i0-1) and x = i1 + ½ (back, i = i1).
        // But the spurious contribution stems from the curl reading
        // E_z[i0,*,*] or E_y[i0,*,*] (which sit at x = i0); the natural
        // sample point for E_inc on the i-face is therefore x = i0·dx.
        // ---------------------------------------------------------------
        if self.i0 >= 1 && self.i1 < grid.nx {
            let x_front = self.i0 as f64 * dx;
            let x_back = self.i1 as f64 * dx;
            // H_y[i,j,k]: (j ∈ [j0, j1], k ∈ [k0, k1]) cross-section.
            // E_inc_z sampled at Ez[i_face, j, k] = (i_face·dx, j·dy, (k+½)·dz).
            let k_hi = self.k1.min(grid.nz.saturating_sub(1));
            let j_hi = self.j1.min(grid.ny);
            for j in self.j0..=j_hi {
                let y = j as f64 * dy;
                for k in self.k0..=k_hi {
                    let z = (k as f64 + 0.5) * dz;
                    let einc_z_f = self.e_inc_component(2, x_front, y, z);
                    let einc_z_b = self.e_inc_component(2, x_back, y, z);
                    grid.hy[(self.i0 - 1, j, k)] -= cx * einc_z_f;
                    grid.hy[(self.i1, j, k)] += cx * einc_z_b;
                }
            }
            // H_z[i,j,k] cross-section: (j ∈ [j0, j1-1], k ∈ [k0, k1+1])
            // — H_z at (i+½, j+½, k); j+½ must lie in [j0, j1]
            // → j ∈ [j0, j1-1]. E_inc_y sampled at Ey[i_face,j,k] =
            // (i_face·dx, (j+½)·dy, k·dz).
            // H_z has shape [nx, ny, nz+1].
            //
            // Phase 2.fdtd.5.3.2: k range is [k0..=k1+1] (not [k0..=k1])
            // — the legacy z-convention places the TF E_y back face at
            // z = (k1+1)·dz (Ey[*, *, k1+1] is TF on the boundary), so
            // the H_z[*, *, k1+1] curl read of E_y[*, *, k1+1] also
            // straddles the i-face when i = i0 - 1 or i = i1. Missing
            // this row contributed materially to the oblique 30°/45°
            // hi-z SF leakage.
            if self.j1 >= 1 {
                let j_hi_hz = self.j1.saturating_sub(1).min(grid.ny.saturating_sub(1));
                let k_hi_hz = (self.k1 + 1).min(grid.nz);
                for j in self.j0..=j_hi_hz {
                    let y = (j as f64 + 0.5) * dy;
                    for k in self.k0..=k_hi_hz {
                        let z = k as f64 * dz;
                        let einc_y_f = self.e_inc_component(1, x_front, y, z);
                        let einc_y_b = self.e_inc_component(1, x_back, y, z);
                        grid.hz[(self.i0 - 1, j, k)] += cx * einc_y_f;
                        grid.hz[(self.i1, j, k)] -= cx * einc_y_b;
                    }
                }
            }
        }

        // ---------------------------------------------------------------
        // j-faces: affect H_x at j0-1, j1 (uses E_inc_z) and
        //          H_z at j0-1, j1 (uses E_inc_x).
        // ---------------------------------------------------------------
        if self.j0 >= 1 && self.j1 < grid.ny {
            let y_front = self.j0 as f64 * dy;
            let y_back = self.j1 as f64 * dy;
            // H_x[i,j,k] cross-section: (i ∈ [i0, i1], k ∈ [k0, k1])
            // H_x at (i, j+½, k+½). E_inc_z at Ez[i, j_face, k] =
            // (i·dx, j_face·dy, (k+½)·dz).
            let i_hi = self.i1.min(grid.nx);
            let k_hi = self.k1.min(grid.nz.saturating_sub(1));
            for i in self.i0..=i_hi {
                let x = i as f64 * dx;
                for k in self.k0..=k_hi {
                    let z = (k as f64 + 0.5) * dz;
                    let einc_z_f = self.e_inc_component(2, x, y_front, z);
                    let einc_z_b = self.e_inc_component(2, x, y_back, z);
                    grid.hx[(i, self.j0 - 1, k)] += cy * einc_z_f;
                    grid.hx[(i, self.j1, k)] -= cy * einc_z_b;
                }
            }
            // H_z[i,j,k] cross-section across j-face: (i ∈ [i0, i1-1],
            // k ∈ [k0, k1+1]). E_inc_x at Ex[i, j_face, k] =
            // ((i+½)·dx, j_face·dy, k·dz).
            // H_z shape [nx, ny, nz+1].
            //
            // Phase 2.fdtd.5.3.2: k range extended to [k0..=k1+1] for
            // the same reason as the i-face H_z block above (TF E_x
            // back face is at z = (k1+1)·dz in the legacy z-convention).
            if self.i1 >= 1 {
                let i_hi_hz = self.i1.saturating_sub(1).min(grid.nx.saturating_sub(1));
                let k_hi_hz = (self.k1 + 1).min(grid.nz);
                for i in self.i0..=i_hi_hz {
                    let x = (i as f64 + 0.5) * dx;
                    for k in self.k0..=k_hi_hz {
                        let z = k as f64 * dz;
                        let einc_x_f = self.e_inc_component(0, x, y_front, z);
                        let einc_x_b = self.e_inc_component(0, x, y_back, z);
                        grid.hz[(i, self.j0 - 1, k)] -= cy * einc_x_f;
                        grid.hz[(i, self.j1, k)] += cy * einc_x_b;
                    }
                }
            }
        }

        // ---------------------------------------------------------------
        // k-faces: the 5.2 z-convention places the back TF boundary at
        // z = (k1+1)·dz (NOT k1·dz; see the front-z vs back-z
        // asymmetry note in `correct_e_kface_plus_x`). So the SF H
        // samples sit at k = k0-1 (front, z = k0-½) and k = k1+1
        // (back, z = (k1+1.5)·dz), and the E_inc sample on the TF
        // side is at the box-boundary plane (z = k0 front, z = (k1+1)
        // back).
        //
        // Affected: H_x at k0-1, k1+1 (uses E_inc_y) and
        //          H_y at k0-1, k1+1 (uses E_inc_x).
        // ---------------------------------------------------------------
        if self.k0 >= 1 && self.k1 + 1 < grid.nz {
            let z_front = self.k0 as f64 * dz;
            let z_back = (self.k1 as f64 + 1.0) * dz;
            // H_x[i,j,k] cross-section across k-face: (i ∈ [i0, i1],
            // j ∈ [j0, j1-1]). E_inc_y at Ey[i, j, k_face] =
            // (i·dx, (j+½)·dy, k_face·dz). H_x shape [nx+1, ny, nz].
            let i_hi = self.i1.min(grid.nx);
            if self.j1 >= 1 {
                let j_hi = self.j1.saturating_sub(1).min(grid.ny.saturating_sub(1));
                for i in self.i0..=i_hi {
                    let x = i as f64 * dx;
                    for j in self.j0..=j_hi {
                        let y = (j as f64 + 0.5) * dy;
                        let einc_y_f = self.e_inc_component(1, x, y, z_front);
                        let einc_y_b = self.e_inc_component(1, x, y, z_back);
                        grid.hx[(i, j, self.k0 - 1)] -= cz * einc_y_f;
                        grid.hx[(i, j, self.k1 + 1)] += cz * einc_y_b;
                    }
                }
            }
            // H_y[i,j,k] cross-section across k-face: (i ∈ [i0, i1-1],
            // j ∈ [j0, j1]). E_inc_x at Ex[i, j, k_face] =
            // ((i+½)·dx, j·dy, k_face·dz). H_y shape [nx, ny+1, nz].
            if self.i1 >= 1 {
                let i_hi_hy = self.i1.saturating_sub(1).min(grid.nx.saturating_sub(1));
                let j_hi_hy = self.j1.min(grid.ny);
                for i in self.i0..=i_hi_hy {
                    let x = (i as f64 + 0.5) * dx;
                    for j in self.j0..=j_hi_hy {
                        let y = j as f64 * dy;
                        let einc_x_f = self.e_inc_component(0, x, y, z_front);
                        let einc_x_b = self.e_inc_component(0, x, y, z_back);
                        grid.hy[(i, j, self.k0 - 1)] += cz * einc_x_f;
                        grid.hy[(i, j, self.k1 + 1)] -= cz * einc_x_b;
                    }
                }
            }
        }
    }

    /// Apply the full 12-stencil TF/SF electric-field correction for
    /// an oblique plane wave. Mirror structure of
    /// [`Self::correct_h_oblique`].
    fn correct_e_oblique(&self, grid: &mut YeeGrid) {
        let dx = grid.dx;
        let dy = grid.dy;
        let dz = grid.dz;
        let dt = self.dt;
        let cx = dt / (EPS0 * dx);
        let cy = dt / (EPS0 * dy);
        let cz = dt / (EPS0 * dz);

        // ---------------------------------------------------------------
        // i-faces: affect E_y at i0, i1 (uses H_inc_z across i-face)
        //          and E_z at i0, i1 (uses H_inc_y across i-face).
        // E_y[i,j,k] at (i, j+½, k); E_z[i,j,k] at (i, j+½, k+½) — wait,
        // E_z[i,j,k] at (i, j, k+½). The straddling H sample on the
        // front side is at i = i0-1 (SF), back is at i = i1 (SF). The
        // H_inc sample position is x = (i0 - ½)·dx for the front and
        // x = (i1 + ½)·dx for the back. But since the aux 1-D grid uses
        // the standard E/H staggering, sampling at (i0-½)·dx is
        // exactly inc_h[pad - 1] for normal incidence; for oblique we
        // project the half-cell-offset x onto k_hat.
        // ---------------------------------------------------------------
        if self.i0 >= 1 && self.i1 < grid.nx {
            // Front i-face sample point for H: x = (i0 - ½)·dx
            // (Hz[i0-1, *, *] sits at x = (i0-1)+½ = i0-½).
            let x_front = (self.i0 as f64 - 0.5) * dx;
            let x_back = (self.i1 as f64 + 0.5) * dx;
            // E_z[i_face, j, k]: (j ∈ [j0, j1], k ∈ [k0, k1])
            // E_z at (i, j, k+½); for samples we use the H_inc_y at
            // (x_face, j·dy, (k+½)·dz).
            let k_hi = self.k1.min(grid.nz.saturating_sub(1));
            let j_hi = self.j1.min(grid.ny);
            for j in self.j0..=j_hi {
                let y = j as f64 * dy;
                for k in self.k0..=k_hi {
                    let z = (k as f64 + 0.5) * dz;
                    let hinc_y_f = self.h_inc_component(1, x_front, y, z);
                    let hinc_y_b = self.h_inc_component(1, x_back, y, z);
                    grid.ez[(self.i0, j, k)] -= cx * hinc_y_f;
                    grid.ez[(self.i1, j, k)] += cx * hinc_y_b;
                }
            }
            // E_y[i_face, j, k]: (j ∈ [j0, j1-1], k ∈ [k0, k1+1])
            // E_y at (i, j+½, k); H_inc_z sample at (x_face, (j+½)·dy, k·dz).
            if self.j1 >= 1 {
                let j_hi_ey = self.j1.saturating_sub(1).min(grid.ny.saturating_sub(1));
                let k_hi_ey = (self.k1 + 1).min(grid.nz);
                for j in self.j0..=j_hi_ey {
                    let y = (j as f64 + 0.5) * dy;
                    for k in self.k0..=k_hi_ey {
                        let z = k as f64 * dz;
                        let hinc_z_f = self.h_inc_component(2, x_front, y, z);
                        let hinc_z_b = self.h_inc_component(2, x_back, y, z);
                        grid.ey[(self.i0, j, k)] += cx * hinc_z_f;
                        grid.ey[(self.i1, j, k)] -= cx * hinc_z_b;
                    }
                }
            }
        }

        // ---------------------------------------------------------------
        // j-faces: affect E_x at j0, j1 (uses H_inc_z across j-face)
        //          and E_z at j0, j1 (uses H_inc_x across j-face).
        // ---------------------------------------------------------------
        if self.j0 >= 1 && self.j1 < grid.ny {
            let y_front = (self.j0 as f64 - 0.5) * dy;
            let y_back = (self.j1 as f64 + 0.5) * dy;
            // E_z[i, j_face, k]: (i ∈ [i0, i1], k ∈ [k0, k1])
            // H_inc_x at (i·dx, y_face, (k+½)·dz).
            let i_hi = self.i1.min(grid.nx);
            let k_hi = self.k1.min(grid.nz.saturating_sub(1));
            for i in self.i0..=i_hi {
                let x = i as f64 * dx;
                for k in self.k0..=k_hi {
                    let z = (k as f64 + 0.5) * dz;
                    let hinc_x_f = self.h_inc_component(0, x, y_front, z);
                    let hinc_x_b = self.h_inc_component(0, x, y_back, z);
                    grid.ez[(i, self.j0, k)] += cy * hinc_x_f;
                    grid.ez[(i, self.j1, k)] -= cy * hinc_x_b;
                }
            }
            // E_x[i, j_face, k]: (i ∈ [i0, i1-1], k ∈ [k0, k1+1])
            // E_x at ((i+½)·dx, j_face·dy, k·dz); H_inc_z sample at
            // ((i+½)·dx, y_face, k·dz). E_x shape [nx, ny+1, nz+1].
            if self.i1 >= 1 {
                let i_hi_ex = self.i1.saturating_sub(1).min(grid.nx.saturating_sub(1));
                let k_hi_ex = (self.k1 + 1).min(grid.nz);
                for i in self.i0..=i_hi_ex {
                    let x = (i as f64 + 0.5) * dx;
                    for k in self.k0..=k_hi_ex {
                        let z = k as f64 * dz;
                        let hinc_z_f = self.h_inc_component(2, x, y_front, z);
                        let hinc_z_b = self.h_inc_component(2, x, y_back, z);
                        grid.ex[(i, self.j0, k)] -= cy * hinc_z_f;
                        grid.ex[(i, self.j1, k)] += cy * hinc_z_b;
                    }
                }
            }
        }

        // ---------------------------------------------------------------
        // k-faces: affect E_x at k0, k1+1 (uses H_inc_y across k-face)
        //          and E_y at k0, k1+1 (uses H_inc_x across k-face).
        // Note the back-face index is k1+1, not k1 — this is the 5.2
        // convention for E components on the top z-face of the TF box.
        // ---------------------------------------------------------------
        if self.k0 >= 1 && self.k1 < grid.nz {
            // k-face z-convention: TF E boundary at z = k0 (front) and
            // z = (k1+1) (back). SF H samples sit at z = (k0-½) (front,
            // H index k = k0-1) and z = (k1+1.5) (back, H index
            // k = k1+1).
            let z_front = (self.k0 as f64 - 0.5) * dz;
            let z_back = (self.k1 as f64 + 1.5) * dz;
            // E_x[i, j, k_face_e]: (i ∈ [i0, i1-1], j ∈ [j0, j1])
            // E_x at ((i+½)·dx, j·dy, k_face_e·dz); H_inc_y sample at
            // ((i+½)·dx, j·dy, z_face).
            if self.i1 >= 1 {
                let i_hi_ex = self.i1.saturating_sub(1).min(grid.nx.saturating_sub(1));
                let j_hi_ex = self.j1.min(grid.ny);
                for i in self.i0..=i_hi_ex {
                    let x = (i as f64 + 0.5) * dx;
                    for j in self.j0..=j_hi_ex {
                        let y = j as f64 * dy;
                        let hinc_y_f = self.h_inc_component(1, x, y, z_front);
                        let hinc_y_b = self.h_inc_component(1, x, y, z_back);
                        grid.ex[(i, j, self.k0)] += cz * hinc_y_f;
                        grid.ex[(i, j, self.k1 + 1)] -= cz * hinc_y_b;
                    }
                }
            }
            // E_y[i, j, k_face_e]: (i ∈ [i0, i1], j ∈ [j0, j1-1])
            // E_y at (i·dx, (j+½)·dy, k_face_e·dz); H_inc_x at
            // (i·dx, (j+½)·dy, z_face). E_y shape [nx+1, ny, nz+1].
            if self.j1 >= 1 {
                let i_hi_ey = self.i1.min(grid.nx);
                let j_hi_ey = self.j1.saturating_sub(1).min(grid.ny.saturating_sub(1));
                for i in self.i0..=i_hi_ey {
                    let x = i as f64 * dx;
                    for j in self.j0..=j_hi_ey {
                        let y = (j as f64 + 0.5) * dy;
                        let hinc_x_f = self.h_inc_component(0, x, y, z_front);
                        let hinc_x_b = self.h_inc_component(0, x, y, z_back);
                        grid.ey[(i, j, self.k0)] -= cz * hinc_x_f;
                        grid.ey[(i, j, self.k1 + 1)] += cz * hinc_x_b;
                    }
                }
            }
        }
    }

    /// Read access to the auxiliary 1-D incident-E grid (mostly for tests).
    pub fn inc_e(&self) -> &[f64] {
        &self.inc_e
    }

    /// Read access to the auxiliary 1-D incident-H grid (mostly for tests).
    pub fn inc_h(&self) -> &[f64] {
        &self.inc_h
    }

    /// Current 1-D step counter.
    pub fn step_count(&self) -> usize {
        self.step
    }
}

// ---- legacy point-source helpers (Phase 2.0) ----

/// Add a Gaussian-time pulse to `E_z(i, j, k)`.
///
/// The injected value is `exp(-((t - t0) / sigma)²)` (a unit-amplitude soft
/// source). The caller controls the time stepping; this function simply
/// *adds* the source contribution to the existing field value.
///
/// # Panics
///
/// Panics if `(i, j, k)` is outside the bounds of `E_z`
/// (shape `[nx+1, ny+1, nz]`).
pub fn gaussian_pulse_ez(
    grid: &mut YeeGrid,
    i: usize,
    j: usize,
    k: usize,
    t: f64,
    t0: f64,
    sigma: f64,
) {
    assert!(
        sigma > 0.0 && sigma.is_finite(),
        "gaussian sigma must be positive and finite"
    );
    let arg = (t - t0) / sigma;
    let amplitude = (-arg * arg).exp();
    grid.ez[(i, j, k)] += amplitude;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Bisection sanity test for [`compute_aux_step`].
    ///
    /// For θ=30°, φ=45°, dx=5 mm, f=3 GHz with a 0.9-Courant dt
    /// (the [`YeeGrid::vacuum`] default), the dispersion-matched
    /// aux step should land within `(0.5·dx, 2.0·dx)` and the
    /// residual `f(ds_aux)` of the dispersion-matching equation
    /// must satisfy `|f| < 1e-12`.
    #[test]
    fn compute_aux_step_root_is_in_bracket_and_residual_is_tight() {
        let theta = 30.0_f64.to_radians();
        let phi = 45.0_f64.to_radians();
        let dx = 5.0e-3_f64;
        let freq = 3.0e9_f64;
        let omega = TAU * freq;
        // dt = 0.9 · dx / (c · √3) matches YeeGrid::vacuum's Courant choice.
        let dt = 0.9 * dx / (C0 * 3.0_f64.sqrt());

        let ds_aux = compute_aux_step(theta, phi, dx, C0, omega, dt);

        // Bracket check.
        assert!(
            ds_aux > 0.5 * dx && ds_aux < 2.0 * dx,
            "ds_aux={ds_aux:.6e} is outside the (0.5·dx, 2.0·dx) primary bracket (dx={dx:.3e})"
        );

        // Residual check: |sin(k_3D · ds/2)/ds − sin(ω·dt/2)/(c·dt)| < 1e-12.
        let k_3d = solve_k_3d(theta, phi, dx, C0, omega, dt);
        let rhs = (omega * dt / 2.0).sin() / (C0 * dt);
        let residual = ((k_3d * ds_aux / 2.0).sin() / ds_aux - rhs).abs();
        assert!(
            residual < 1.0e-12,
            "dispersion-match residual {residual:.3e} exceeds 1e-12 tolerance"
        );

        // For this specific case, ds_aux should be ~0.770·dx (independent of
        // dt as long as it's CFL-stable). Use a loose 5% window.
        let ratio = ds_aux / dx;
        assert!(
            (0.7..0.85).contains(&ratio),
            "ds_aux/dx={ratio:.4} far from the expected ~0.77 for θ=30°, φ=45°"
        );
    }

    /// On-axis incidence (θ=π/2, φ=0) makes the 1-D and 3-D Yee
    /// dispersion relations coincide exactly: `k_3D · dx = ω·dt · (dx/(c·dt))`
    /// once the only nonzero `k̂_α` lands on a single coordinate axis,
    /// and the dispersion-matched `ds_aux` collapses to `dx`. Tight
    /// tolerance.
    #[test]
    fn compute_aux_step_collapses_to_dx_on_axis() {
        let theta = PI / 2.0; // +xy plane
        let phi = 0.0; // +x
        let dx = 5.0e-3_f64;
        let freq = 3.0e9_f64;
        let omega = TAU * freq;
        let dt = 0.9 * dx / (C0 * 3.0_f64.sqrt());

        let ds_aux = compute_aux_step(theta, phi, dx, C0, omega, dt);
        assert!(
            (ds_aux - dx).abs() < 1.0e-9 * dx,
            "on-axis ds_aux ({ds_aux:.6e}) should collapse to dx ({dx:.6e})"
        );
    }
}
