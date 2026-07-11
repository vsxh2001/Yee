//! Uniform-vacuum FDTD problem description (E.0 scope).

use yee_core::units::C0;

/// Description of a uniform, lossless FDTD problem on a Yee grid.
///
/// This is the E.0 walking-skeleton contract (ADR-0175): scalar `eps_r` /
/// `mu_r`, σ = 0, PEC outer box. The staggered component shapes follow
/// `yee_fdtd::grid::YeeGrid` exactly so parity is testable index-for-index.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FdtdSpec {
    /// Number of primary cells along x.
    pub nx: usize,
    /// Number of primary cells along y.
    pub ny: usize,
    /// Number of primary cells along z.
    pub nz: usize,
    /// Cell size along x (meters).
    pub dx: f64,
    /// Cell size along y (meters).
    pub dy: f64,
    /// Cell size along z (meters).
    pub dz: f64,
    /// Time step (seconds). Must be ≤ [`FdtdSpec::courant_limit`].
    pub dt: f64,
    /// Relative permittivity (uniform; vacuum = 1.0).
    pub eps_r: f64,
    /// Relative permeability (uniform; vacuum = 1.0).
    pub mu_r: f64,
}

impl FdtdSpec {
    /// Build a vacuum cubic-cell spec with `dt` at 0.9× the Courant limit,
    /// mirroring `YeeGrid::vacuum`.
    ///
    /// # Panics
    ///
    /// Panics if any dimension is zero or `dx` is non-positive/non-finite.
    pub fn vacuum(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        assert!(nx > 0 && ny > 0 && nz > 0, "grid dimensions must be > 0");
        assert!(
            dx.is_finite() && dx > 0.0,
            "cell size must be positive and finite"
        );
        let mut spec = Self {
            nx,
            ny,
            nz,
            dx,
            dy: dx,
            dz: dx,
            dt: 0.0,
            eps_r: 1.0,
            mu_r: 1.0,
        };
        spec.dt = 0.9 * spec.courant_limit();
        spec
    }

    /// CFL stability limit `1 / (c₀ √(1/dx² + 1/dy² + 1/dz²))`.
    pub fn courant_limit(&self) -> f64 {
        1.0 / (C0
            * (1.0 / (self.dx * self.dx) + 1.0 / (self.dy * self.dy) + 1.0 / (self.dz * self.dz))
                .sqrt())
    }

    /// Shape of the `E_x` array: `[nx, ny+1, nz+1]`.
    pub fn ex_dims(&self) -> (usize, usize, usize) {
        (self.nx, self.ny + 1, self.nz + 1)
    }

    /// Shape of the `E_y` array: `[nx+1, ny, nz+1]`.
    pub fn ey_dims(&self) -> (usize, usize, usize) {
        (self.nx + 1, self.ny, self.nz + 1)
    }

    /// Shape of the `E_z` array: `[nx+1, ny+1, nz]`.
    pub fn ez_dims(&self) -> (usize, usize, usize) {
        (self.nx + 1, self.ny + 1, self.nz)
    }

    /// Shape of the `H_x` array: `[nx+1, ny, nz]`.
    pub fn hx_dims(&self) -> (usize, usize, usize) {
        (self.nx + 1, self.ny, self.nz)
    }

    /// Shape of the `H_y` array: `[nx, ny+1, nz]`.
    pub fn hy_dims(&self) -> (usize, usize, usize) {
        (self.nx, self.ny + 1, self.nz)
    }

    /// Shape of the `H_z` array: `[nx, ny, nz+1]`.
    pub fn hz_dims(&self) -> (usize, usize, usize) {
        (self.nx, self.ny, self.nz + 1)
    }
}

/// Per-axis nonuniform **primal cell widths** (FS.0b.0, ADR-0208).
///
/// `dx`/`dy`/`dz` hold one width (metres) per primary cell, so their lengths
/// must equal `nx`/`ny`/`nz`. Constant arrays reproduce the uniform kernel
/// **bit-exactly** (gate `compute-018`): the kernel divides by these values
/// (and by the derived dual spacings, which are bit-equal to the primal for
/// constant arrays), never by precomputed inverses.
#[derive(Debug, Clone, PartialEq)]
pub struct GradedSpacings {
    /// Primal cell widths along x (length `nx`, metres).
    pub dx: Vec<f64>,
    /// Primal cell widths along y (length `ny`, metres).
    pub dy: Vec<f64>,
    /// Primal cell widths along z (length `nz`, metres).
    pub dz: Vec<f64>,
}

impl GradedSpacings {
    /// Check array lengths against `spec` and every width for positivity /
    /// finiteness. Returns a human-readable error (job specs arrive over
    /// untrusted transports; the engine forwards this as an error event).
    pub fn validate(&self, spec: &FdtdSpec) -> Result<(), String> {
        for (axis, arr, n) in [
            ("dx", &self.dx, spec.nx),
            ("dy", &self.dy, spec.ny),
            ("dz", &self.dz, spec.nz),
        ] {
            if arr.len() != n {
                return Err(format!(
                    "spacings.{axis} has {} entries, expected {n} (one primal cell width per cell)",
                    arr.len()
                ));
            }
            if let Some((i, &d)) = arr
                .iter()
                .enumerate()
                .find(|(_, d)| !(d.is_finite() && **d > 0.0))
            {
                return Err(format!(
                    "spacings.{axis}[{i}] = {d} must be positive and finite"
                ));
            }
        }
        Ok(())
    }

    /// FS.0b.0 scope rule: spacing must be **uniform within the CPML layers**
    /// of every absorbing face (`faces` is
    /// `[[x−, x+], [y−, y+], [z−, z+]]`). The Roden–Gedney profile grading
    /// assumes one cell size per layer; mesh rules never grade inside
    /// absorbers. Call after [`GradedSpacings::validate`].
    pub fn validate_cpml_layers(&self, npml: usize, faces: [[bool; 2]; 3]) -> Result<(), String> {
        for (axis, arr, f) in [
            ("dx", &self.dx, faces[0]),
            ("dy", &self.dy, faces[1]),
            ("dz", &self.dz, faces[2]),
        ] {
            let layer = npml.min(arr.len());
            if layer == 0 {
                continue;
            }
            if f[0] && arr[..layer].iter().any(|d| *d != arr[0]) {
                return Err(format!(
                    "spacings.{axis} is graded inside the {axis}-min CPML layer \
                     (the first {layer} cells must share one width; FS.0b.0 scope)"
                ));
            }
            let hi = &arr[arr.len() - layer..];
            if f[1] && hi.iter().any(|d| *d != hi[0]) {
                return Err(format!(
                    "spacings.{axis} is graded inside the {axis}-max CPML layer \
                     (the last {layer} cells must share one width; FS.0b.0 scope)"
                ));
            }
        }
        Ok(())
    }

    /// CFL stability limit from the **minimum** spacing per axis — the same
    /// expression shape as [`FdtdSpec::courant_limit`], so constant arrays
    /// yield a bit-identical limit.
    ///
    /// # Panics
    ///
    /// Panics if any axis array is empty (rejected by
    /// [`GradedSpacings::validate`]).
    pub fn courant_limit(&self) -> f64 {
        let min = |arr: &[f64]| arr.iter().copied().fold(f64::INFINITY, f64::min);
        let (dx, dy, dz) = (min(&self.dx), min(&self.dy), min(&self.dz));
        assert!(
            dx.is_finite() && dy.is_finite() && dz.is_finite(),
            "GradedSpacings::courant_limit on an empty axis"
        );
        1.0 / (C0 * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)).sqrt())
    }
}

/// One axis's kernel-ready spacings: `primal[i]` is the width of primary
/// cell `i` (length `n`); `dual[i]` is the distance between the H samples
/// straddling node `i` — `(primal[i−1] + primal[i]) / 2` in the interior,
/// the single adjacent primal at the domain edges (length `n + 1`).
#[derive(Debug, Clone)]
pub(crate) struct AxisSpacings {
    pub(crate) primal: Vec<f64>,
    pub(crate) dual: Vec<f64>,
}

impl AxisSpacings {
    /// Uniform fill: every primal and dual entry is bit-equal to `d`.
    fn uniform(n: usize, d: f64) -> Self {
        Self {
            primal: vec![d; n],
            dual: vec![d; n + 1],
        }
    }

    /// Build the dual array from primal widths. For a constant array this
    /// reduces to the uniform fill bit-exactly (`(d + d) / 2 == d` in
    /// IEEE-754: exact ×2 then exact halving).
    fn from_primal(primal: Vec<f64>) -> Self {
        let n = primal.len();
        let mut dual = Vec::with_capacity(n + 1);
        dual.push(primal[0]);
        for i in 1..n {
            dual.push((primal[i - 1] + primal[i]) / 2.0);
        }
        dual.push(primal[n - 1]);
        Self { primal, dual }
    }
}

/// Kernel-ready spacings for all three axes. The kernel always divides by
/// these entries; the uniform constructors fill them with the scalar
/// `spec.dx/dy/dz`, so the uniform path is bit-exact by construction
/// (division is a pure function of its operand bit patterns).
#[derive(Debug, Clone)]
pub(crate) struct SpacingArrays {
    pub(crate) x: AxisSpacings,
    pub(crate) y: AxisSpacings,
    pub(crate) z: AxisSpacings,
}

impl SpacingArrays {
    /// Uniform spacings from the scalar spec (today's path).
    pub(crate) fn uniform(spec: &FdtdSpec) -> Self {
        Self {
            x: AxisSpacings::uniform(spec.nx, spec.dx),
            y: AxisSpacings::uniform(spec.ny, spec.dy),
            z: AxisSpacings::uniform(spec.nz, spec.dz),
        }
    }

    /// Graded spacings from validated primal widths.
    pub(crate) fn graded(graded: &GradedSpacings) -> Self {
        Self {
            x: AxisSpacings::from_primal(graded.dx.clone()),
            y: AxisSpacings::from_primal(graded.dy.clone()),
            z: AxisSpacings::from_primal(graded.dz.clone()),
        }
    }
}

/// Flat row-major index into an array of shape `dims`, matching `ndarray`'s
/// default (C-order) layout: `(i * dim_j + j) * dim_k + k`.
#[inline]
pub(crate) fn idx3(dims: (usize, usize, usize), i: usize, j: usize, k: usize) -> usize {
    debug_assert!(i < dims.0 && j < dims.1 && k < dims.2);
    (i * dims.1 + j) * dims.2 + k
}

/// Element count of an array of shape `dims`.
#[inline]
pub(crate) fn len3(dims: (usize, usize, usize)) -> usize {
    dims.0 * dims.1 * dims.2
}

#[cfg(test)]
mod tests {
    use super::*;

    fn constant(spec: &FdtdSpec) -> GradedSpacings {
        GradedSpacings {
            dx: vec![spec.dx; spec.nx],
            dy: vec![spec.dy; spec.ny],
            dz: vec![spec.dz; spec.nz],
        }
    }

    #[test]
    fn dual_spacings_interior_and_edges() {
        // Known graded axis: primal [1, 2, 4] → dual [1, 1.5, 3, 4].
        let ax = AxisSpacings::from_primal(vec![1.0, 2.0, 4.0]);
        assert_eq!(ax.primal, vec![1.0, 2.0, 4.0]);
        assert_eq!(ax.dual, vec![1.0, 1.5, 3.0, 4.0]);
    }

    #[test]
    fn constant_arrays_reduce_to_uniform_bit_exactly() {
        let spec = FdtdSpec::vacuum(5, 7, 3, 0.37e-3);
        let g = constant(&spec);
        let graded = SpacingArrays::graded(&g);
        let uniform = SpacingArrays::uniform(&spec);
        for (a, b) in [
            (&graded.x, &uniform.x),
            (&graded.y, &uniform.y),
            (&graded.z, &uniform.z),
        ] {
            // Exact f64 equality — the bit-exact-on-uniform foundation.
            assert!(a.primal.iter().zip(&b.primal).all(|(p, q)| p == q));
            assert!(a.dual.iter().zip(&b.dual).all(|(p, q)| p == q));
            assert_eq!(a.dual.len(), a.primal.len() + 1);
        }
        // Courant limit from constant arrays is bit-identical to the scalar.
        assert!(g.courant_limit() == spec.courant_limit());
    }

    #[test]
    fn courant_limit_uses_minimum_spacing() {
        let spec = FdtdSpec::vacuum(4, 4, 4, 1.0e-3);
        let mut g = constant(&spec);
        g.dx[2] = 0.25e-3; // one fine cell tightens the limit
        let fine = FdtdSpec {
            dx: 0.25e-3,
            ..spec
        };
        assert!(g.courant_limit() == fine.courant_limit());
    }

    #[test]
    fn validate_rejects_bad_lengths_and_values() {
        let spec = FdtdSpec::vacuum(4, 4, 4, 1.0e-3);
        let mut g = constant(&spec);
        g.dy.pop();
        let err = g.validate(&spec).unwrap_err();
        assert!(err.contains("spacings.dy"), "{err}");
        assert!(err.contains("expected 4"), "{err}");

        let mut g = constant(&spec);
        g.dz[1] = 0.0;
        let err = g.validate(&spec).unwrap_err();
        assert!(err.contains("spacings.dz[1]"), "{err}");

        let mut g = constant(&spec);
        g.dx[3] = f64::NAN;
        assert!(g.validate(&spec).is_err());

        assert!(constant(&spec).validate(&spec).is_ok());
    }

    #[test]
    fn cpml_layers_must_be_uniform_on_absorbing_faces() {
        let spec = FdtdSpec::vacuum(12, 12, 12, 1.0e-3);
        let g = constant(&spec);
        assert!(g.validate_cpml_layers(4, [[true; 2]; 3]).is_ok());

        // Grading inside the x-min layer is rejected …
        let mut bad = constant(&spec);
        bad.dx[1] = 0.5e-3;
        let err = bad.validate_cpml_layers(4, [[true; 2]; 3]).unwrap_err();
        assert!(err.contains("x-min CPML layer"), "{err}");
        // … unless that face is not absorbing.
        assert!(
            bad.validate_cpml_layers(4, [[false, true], [true; 2], [true; 2]])
                .is_ok()
        );

        // Same for the max side.
        let mut bad = constant(&spec);
        bad.dz[11] = 0.5e-3;
        let err = bad.validate_cpml_layers(4, [[true; 2]; 3]).unwrap_err();
        assert!(err.contains("z-max CPML layer"), "{err}");

        // Grading in the interior is fine.
        let mut ok = constant(&spec);
        ok.dy[6] = 0.5e-3;
        assert!(ok.validate_cpml_layers(4, [[true; 2]; 3]).is_ok());
    }
}
