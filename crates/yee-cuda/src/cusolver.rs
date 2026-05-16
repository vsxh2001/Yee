//! Dense complex-double LU factorization + triangular solve via cuSOLVER.
//!
//! Phase 1.5 wires [`cusolverDnZgetrf`] + [`cusolverDnZgetrs`] behind a
//! safe RAII type, [`DenseLuComplex`]. Without `--features cuda`, the type
//! exists but every constructor / method returns [`Error::NotEnabled`].
//!
//! Matrices are passed in **column-major** order, matching LAPACK and the
//! cuSOLVER C ABI:
//!
//! ```text
//! A[i + j*n]  // row i, column j, in a length-(n*n) Vec<Complex64>
//! ```
//!
//! This commit lands the public API + the no-feature stub. The next
//! commit wires up the actual cudarc/cuSOLVER calls behind the `cuda`
//! feature.
//!
//! [`cusolverDnZgetrf`]: https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdngetrf
//! [`cusolverDnZgetrs`]: https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdngetrs

use num_complex::Complex64;

use crate::{Error, Result};

/// Pre-factored dense complex-double LU.
///
/// Holds one instance per matrix you need to solve repeatedly. Each
/// instance owns its device-side LU factors and pivot vector so
/// [`Self::solve`] can be called many times against new right-hand sides
/// without re-factoring.
///
/// # Feature gating
///
/// All constructors return [`Error::NotEnabled`] unless the crate is
/// compiled with `--features cuda`. The struct itself still exists in
/// both configurations so downstream code can name the type
/// unconditionally.
pub struct DenseLuComplex {
    // Zero-sized marker keeps the struct nameable in the no-feature
    // build. The `cuda` build will replace this with the real device
    // state in the next commit.
    _marker: core::marker::PhantomData<Complex64>,
}

impl DenseLuComplex {
    /// Factorize an `n`×`n` complex-double matrix `a` (column-major).
    ///
    /// # Errors
    ///
    /// - [`Error::NotEnabled`] if compiled without `--features cuda`.
    /// - (Phase 1.5b) [`Error::Driver`] if any underlying CUDA / cuSOLVER
    ///   call fails, if `a.len() != n*n`, `n == 0`, or cuSOLVER reports a
    ///   singular `U` factor (`info > 0`).
    pub fn factorize(a: &[Complex64], n: usize) -> Result<Self> {
        let _ = (a, n);
        Err(Error::NotEnabled)
    }

    /// Solve `A x = b` using the precomputed LU.
    ///
    /// `b` is column-major of logical shape `(n, nrhs)`. Returns the
    /// `n`-by-`nrhs` solution `x` in the same layout.
    ///
    /// # Errors
    ///
    /// - [`Error::NotEnabled`] if compiled without `--features cuda`.
    /// - (Phase 1.5b) [`Error::Driver`] if any underlying CUDA / cuSOLVER
    ///   call fails, `b.len() != n*nrhs`, or `nrhs == 0`.
    pub fn solve(&self, b: &[Complex64], nrhs: usize) -> Result<Vec<Complex64>> {
        let _ = (b, nrhs);
        Err(Error::NotEnabled)
    }
}

#[cfg(test)]
mod tests {
    #[cfg(not(feature = "cuda"))]
    use super::*;
    #[cfg(not(feature = "cuda"))]
    use num_complex::Complex64;

    #[cfg(not(feature = "cuda"))]
    #[test]
    fn factorize_without_feature_is_not_enabled() {
        let a = vec![Complex64::new(1.0, 0.0)];
        match DenseLuComplex::factorize(&a, 1) {
            Err(Error::NotEnabled) => {}
            Err(other) => panic!("expected NotEnabled, got {other:?}"),
            Ok(_) => panic!("expected NotEnabled, got Ok(DenseLuComplex)"),
        }
    }
}
