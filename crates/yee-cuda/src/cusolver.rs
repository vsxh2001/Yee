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
//! [`cusolverDnZgetrf`]: https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdngetrf
//! [`cusolverDnZgetrs`]: https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdngetrs

use num_complex::Complex64;

#[cfg(not(feature = "cuda"))]
use crate::Error;
use crate::Result;

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
    #[cfg(feature = "cuda")]
    inner: imp::Inner,
    // Without `cuda`, keep the struct non-empty so future fields can be
    // added without a layout warning. PhantomData is zero-sized.
    #[cfg(not(feature = "cuda"))]
    _marker: core::marker::PhantomData<Complex64>,
}

impl DenseLuComplex {
    /// Factorize an `n`×`n` complex-double matrix `a` (column-major).
    ///
    /// # Errors
    ///
    /// - [`Error::NotEnabled`] if compiled without `--features cuda`.
    /// - [`Error::Driver`] if any underlying CUDA / cuSOLVER call fails,
    ///   if `a.len() != n*n`, `n == 0`, or cuSOLVER reports a singular
    ///   `U` factor (`info > 0`).
    pub fn factorize(a: &[Complex64], n: usize) -> Result<Self> {
        #[cfg(not(feature = "cuda"))]
        {
            let _ = (a, n);
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "cuda")]
        {
            imp::factorize(a, n).map(|inner| Self { inner })
        }
    }

    /// Solve `A x = b` using the precomputed LU.
    ///
    /// `b` is column-major of logical shape `(n, nrhs)`. Returns the
    /// `n`-by-`nrhs` solution `x` in the same layout.
    ///
    /// # Errors
    ///
    /// - [`Error::NotEnabled`] if compiled without `--features cuda`.
    /// - [`Error::Driver`] if any underlying CUDA / cuSOLVER call fails,
    ///   `b.len() != n*nrhs`, or `nrhs == 0`.
    pub fn solve(&self, b: &[Complex64], nrhs: usize) -> Result<Vec<Complex64>> {
        #[cfg(not(feature = "cuda"))]
        {
            let _ = (b, nrhs);
            Err(Error::NotEnabled)
        }
        #[cfg(feature = "cuda")]
        {
            imp::solve(&self.inner, b, nrhs)
        }
    }
}

#[cfg(feature = "cuda")]
#[allow(unsafe_code)]
mod imp {
    //! Cudarc-backed implementation of [`super::DenseLuComplex`].
    //!
    //! The whole module opts back into `unsafe_code` because cuSOLVER's
    //! `sys::*` functions are `unsafe fn`. Each `unsafe` block carries a
    //! local SAFETY comment proving the invariants it relies on.
    use std::sync::Arc;

    use cudarc::cusolver::safe::DnHandle;
    use cudarc::cusolver::sys;
    use cudarc::driver::{CudaContext, CudaSlice, CudaStream, DevicePtr, DevicePtrMut};
    use num_complex::Complex64;

    use crate::{Error, Result};

    /// Owned device-side state for one factored matrix.
    pub(super) struct Inner {
        pub(super) n: usize,
        /// LU factors, stored as `2*n*n` interleaved `f64` (re, im).
        pub(super) lu: CudaSlice<f64>,
        /// Pivot vector (size `n`, 1-indexed in LAPACK convention).
        pub(super) pivots: CudaSlice<i32>,
        /// cuSOLVER dense handle (RAII; destroys on drop).
        pub(super) handle: DnHandle,
        /// Stream the handle was bound to. Kept alive alongside the
        /// handle so `set_stream` does not dangle.
        pub(super) stream: Arc<CudaStream>,
    }

    fn drv<E: std::fmt::Display>(e: E) -> Error {
        Error::Driver(format!("{e}"))
    }

    fn solver<E: std::fmt::Display>(e: E) -> Error {
        Error::Driver(format!("cusolver: {e}"))
    }

    /// Cast `&[Complex64]` -> `&[f64]` of doubled length.
    fn complex_as_f64(src: &[Complex64]) -> &[f64] {
        // SAFETY: `num_complex::Complex<f64>` is `#[repr(C)] { re: f64,
        // im: f64 }` per upstream docs, so a `[Complex<f64>; n]` is
        // bit-identical to `[f64; 2*n]` with the same alignment.
        unsafe { std::slice::from_raw_parts(src.as_ptr().cast::<f64>(), src.len() * 2) }
    }

    /// Cast a `Vec<f64>` of even length into the corresponding
    /// `Vec<Complex64>` without reallocating.
    fn f64_vec_to_complex(mut v: Vec<f64>) -> Vec<Complex64> {
        debug_assert!(v.len().is_multiple_of(2));
        debug_assert!(v.capacity().is_multiple_of(2));
        let len = v.len() / 2;
        let cap = v.capacity() / 2;
        let ptr = v.as_mut_ptr().cast::<Complex64>();
        core::mem::forget(v);
        // SAFETY: `Complex<f64>` is `#[repr(C)] { re: f64, im: f64 }`, so
        // a `[f64; 2*n]` allocation owned by a `Vec<f64>` is layout- and
        // alignment-compatible with `Vec<Complex<f64>>` of half the
        // length/capacity. We forgot the original Vec so there is no
        // double-free. f64 and Complex<f64> share alignment 8.
        unsafe { Vec::from_raw_parts(ptr, len, cap) }
    }

    fn to_cint(v: usize, what: &str) -> Result<core::ffi::c_int> {
        core::ffi::c_int::try_from(v)
            .map_err(|_| Error::Driver(format!("{what}={v} does not fit in c_int")))
    }

    pub(super) fn factorize(a: &[Complex64], n: usize) -> Result<Inner> {
        if n == 0 {
            return Err(Error::Driver("n must be > 0".into()));
        }
        if a.len() != n * n {
            return Err(Error::Driver(format!(
                "factorize: a.len()={} but n*n={}",
                a.len(),
                n * n
            )));
        }
        let n_i = to_cint(n, "n")?;

        let ctx = CudaContext::new(0).map_err(drv)?;
        let stream = ctx.default_stream();
        let handle = DnHandle::new(stream.clone()).map_err(solver)?;

        let mut lu = stream.clone_htod(complex_as_f64(a)).map_err(drv)?;
        let mut pivots = stream.alloc_zeros::<i32>(n).map_err(drv)?;
        let mut info = stream.alloc_zeros::<i32>(1).map_err(drv)?;

        // 1) Workspace query (in cuDoubleComplex elements).
        let mut lwork: core::ffi::c_int = 0;
        {
            let (lu_ptr, _g) = lu.device_ptr_mut(&stream);
            // SAFETY: lu_ptr is a valid device pointer to 2*n*n f64 ==
            // n*n cuDoubleComplex; handle and lda=n_i match the matrix.
            unsafe {
                sys::cusolverDnZgetrf_bufferSize(
                    handle.cu(),
                    n_i,
                    n_i,
                    lu_ptr as *mut sys::cuDoubleComplex,
                    n_i,
                    &mut lwork as *mut core::ffi::c_int,
                )
                .result()
                .map_err(solver)?;
            }
        }
        if lwork < 0 {
            return Err(Error::Driver(format!(
                "cusolverDnZgetrf_bufferSize returned negative lwork={lwork}"
            )));
        }

        // 2) Allocate workspace (lwork elements of cuDoubleComplex ==
        //    2*lwork f64).
        let lwork_usize = lwork as usize;
        let mut work = stream.alloc_zeros::<f64>(2 * lwork_usize).map_err(drv)?;

        // 3) Factor.
        {
            let (lu_ptr, _g_lu) = lu.device_ptr_mut(&stream);
            let (work_ptr, _g_w) = work.device_ptr_mut(&stream);
            let (piv_ptr, _g_p) = pivots.device_ptr_mut(&stream);
            let (info_ptr, _g_i) = info.device_ptr_mut(&stream);
            // SAFETY: all pointers reference correctly-sized device
            // buffers (A: n*n, Workspace: lwork cuDoubleComplex; pivots:
            // n i32; info: 1 i32). The handle was created above.
            unsafe {
                sys::cusolverDnZgetrf(
                    handle.cu(),
                    n_i,
                    n_i,
                    lu_ptr as *mut sys::cuDoubleComplex,
                    n_i,
                    work_ptr as *mut sys::cuDoubleComplex,
                    piv_ptr as *mut core::ffi::c_int,
                    info_ptr as *mut core::ffi::c_int,
                )
                .result()
                .map_err(solver)?;
            }
        }

        stream.synchronize().map_err(drv)?;
        let info_host = stream.clone_dtoh(&info).map_err(drv)?;
        if info_host[0] != 0 {
            return Err(Error::Driver(format!(
                "cusolverDnZgetrf returned info={} (negative=bad arg, positive=singular U)",
                info_host[0]
            )));
        }
        drop(work);

        Ok(Inner {
            n,
            lu,
            pivots,
            handle,
            stream,
        })
    }

    /// One-shot Zgetrf that returns the LU factors and pivots to the host.
    ///
    /// Used by [`crate::backend::Backend::cusolver_zgetrf`]; not part of
    /// the [`super::DenseLuComplex`] flow (which keeps everything on the
    /// device for repeat solves).
    pub(crate) fn zgetrf_host(
        a: &[Complex64],
        n: usize,
    ) -> Result<(Vec<Complex64>, Vec<i32>)> {
        let inner = factorize(a, n)?;
        let lu_f64 = inner.stream.clone_dtoh(&inner.lu).map_err(drv)?;
        let pivots = inner.stream.clone_dtoh(&inner.pivots).map_err(drv)?;
        Ok((f64_vec_to_complex(lu_f64), pivots))
    }

    /// One-shot Zgetrs that uploads `lu`+`pivots` and solves for `x`.
    ///
    /// Used by [`crate::backend::Backend::cusolver_zgetrs`].
    pub(crate) fn zgetrs_host(
        lu: &[Complex64],
        pivots: &[i32],
        b: &[Complex64],
        n: usize,
        nrhs: usize,
    ) -> Result<Vec<Complex64>> {
        if n == 0 {
            return Err(Error::Driver("n must be > 0".into()));
        }
        if lu.len() != n * n {
            return Err(Error::Driver(format!(
                "zgetrs_host: lu.len()={} but n*n={}",
                lu.len(),
                n * n
            )));
        }
        if pivots.len() != n {
            return Err(Error::Driver(format!(
                "zgetrs_host: pivots.len()={} but n={}",
                pivots.len(),
                n
            )));
        }

        let ctx = CudaContext::new(0).map_err(drv)?;
        let stream = ctx.default_stream();
        let handle = DnHandle::new(stream.clone()).map_err(solver)?;

        let lu_dev = stream.clone_htod(complex_as_f64(lu)).map_err(drv)?;
        let piv_dev = stream.clone_htod(pivots).map_err(drv)?;
        let inner = Inner {
            n,
            lu: lu_dev,
            pivots: piv_dev,
            handle,
            stream,
        };
        solve(&inner, b, nrhs)
    }

    pub(super) fn solve(inner: &Inner, b: &[Complex64], nrhs: usize) -> Result<Vec<Complex64>> {
        if nrhs == 0 {
            return Err(Error::Driver("nrhs must be > 0".into()));
        }
        let n = inner.n;
        if b.len() != n * nrhs {
            return Err(Error::Driver(format!(
                "solve: b.len()={} but n*nrhs={}",
                b.len(),
                n * nrhs
            )));
        }
        let n_i = to_cint(n, "n")?;
        let nrhs_i = to_cint(nrhs, "nrhs")?;
        let stream = &inner.stream;

        let mut bx = stream.clone_htod(complex_as_f64(b)).map_err(drv)?;
        let mut info = stream.alloc_zeros::<i32>(1).map_err(drv)?;

        {
            let (lu_ptr, _g_lu) = inner.lu.device_ptr(stream);
            let (piv_ptr, _g_p) = inner.pivots.device_ptr(stream);
            let (bx_ptr, _g_b) = bx.device_ptr_mut(stream);
            let (info_ptr, _g_i) = info.device_ptr_mut(stream);
            // SAFETY: handle from inner is alive; LU and pivots come from
            // a successful Zgetrf on the same n; B is 2*n*nrhs f64 ==
            // n*nrhs cuDoubleComplex; trans=N (no transpose); leading
            // dimensions match the row count.
            unsafe {
                sys::cusolverDnZgetrs(
                    inner.handle.cu(),
                    sys::cublasOperation_t::CUBLAS_OP_N,
                    n_i,
                    nrhs_i,
                    lu_ptr as *const sys::cuDoubleComplex,
                    n_i,
                    piv_ptr as *const core::ffi::c_int,
                    bx_ptr as *mut sys::cuDoubleComplex,
                    n_i,
                    info_ptr as *mut core::ffi::c_int,
                )
                .result()
                .map_err(solver)?;
            }
        }

        stream.synchronize().map_err(drv)?;
        let info_host = stream.clone_dtoh(&info).map_err(drv)?;
        if info_host[0] != 0 {
            return Err(Error::Driver(format!(
                "cusolverDnZgetrs returned info={}",
                info_host[0]
            )));
        }

        let x_f64 = stream.clone_dtoh(&bx).map_err(drv)?;
        Ok(f64_vec_to_complex(x_f64))
    }
}

/// Backend-trait entry point: dense complex-double LU factorization,
/// returning host-side LU factors and pivots. See
/// [`crate::backend::Backend::cusolver_zgetrf`].
#[cfg(feature = "cuda")]
pub(crate) fn zgetrf_host(a: &[Complex64], n: usize) -> Result<(Vec<Complex64>, Vec<i32>)> {
    imp::zgetrf_host(a, n)
}

/// Backend-trait entry point: triangular solve given host-side LU
/// factors and pivots. See [`crate::backend::Backend::cusolver_zgetrs`].
#[cfg(feature = "cuda")]
pub(crate) fn zgetrs_host(
    lu: &[Complex64],
    pivots: &[i32],
    b: &[Complex64],
    n: usize,
    nrhs: usize,
) -> Result<Vec<Complex64>> {
    imp::zgetrs_host(lu, pivots, b, n, nrhs)
}

#[cfg(test)]
mod tests {
    #[cfg(not(feature = "cuda"))]
    use super::*;
    #[cfg(feature = "cuda")]
    use super::DenseLuComplex;
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

    /// Smoke test against the identity. `#[ignore]`-d because it needs a
    /// reachable CUDA driver; run with `cargo test -- --ignored` on a
    /// GPU host.
    #[cfg(feature = "cuda")]
    #[test]
    #[ignore]
    fn factorize_and_solve_2x2_identity() {
        // A = I_2 (column-major), b = [2+3j, 4+5j]. Expect x == b.
        let a = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        ];
        let b = vec![Complex64::new(2.0, 3.0), Complex64::new(4.0, 5.0)];
        let lu = DenseLuComplex::factorize(&a, 2).expect("factorize");
        let x = lu.solve(&b, 1).expect("solve");
        assert!(
            (x[0] - b[0]).norm() < 1e-12,
            "x[0]={:?} vs b[0]={:?}",
            x[0],
            b[0]
        );
        assert!(
            (x[1] - b[1]).norm() < 1e-12,
            "x[1]={:?} vs b[1]={:?}",
            x[1],
            b[1]
        );
    }

    /// 64×64 random Hermitian positive-definite solve.
    ///
    /// Build `A = M^H M + I` (HPD ⇒ non-singular), draw a random `b`,
    /// solve `A x = b`, and verify `||A x − b|| / ||b|| < 1e-10`.
    /// `#[ignore]`-d for the same reason as the 2×2 case.
    #[cfg(feature = "cuda")]
    #[test]
    #[ignore]
    fn factorize_and_solve_random_64() {
        use rand::{Rng, SeedableRng, rngs::StdRng};

        const N: usize = 64;
        let mut rng = StdRng::seed_from_u64(0xfeed_f00d);

        // Random M (column-major n×n).
        let mut m = vec![Complex64::new(0.0, 0.0); N * N];
        for v in &mut m {
            *v = Complex64::new(rng.random_range(-1.0..1.0), rng.random_range(-1.0..1.0));
        }

        // A = M^H * M + I (column-major).
        let mut a = vec![Complex64::new(0.0, 0.0); N * N];
        for j in 0..N {
            for i in 0..N {
                // A[i,j] = sum_k conj(M[k,i]) * M[k,j]
                let mut acc = Complex64::new(0.0, 0.0);
                for k in 0..N {
                    acc += m[k + i * N].conj() * m[k + j * N];
                }
                a[i + j * N] = acc;
                if i == j {
                    a[i + j * N] += Complex64::new(1.0, 0.0);
                }
            }
        }

        // Random b.
        let mut b = vec![Complex64::new(0.0, 0.0); N];
        for v in &mut b {
            *v = Complex64::new(rng.random_range(-1.0..1.0), rng.random_range(-1.0..1.0));
        }

        let lu = DenseLuComplex::factorize(&a, N).expect("factorize");
        let x = lu.solve(&b, 1).expect("solve");

        // residual = A x − b (host).
        let mut residual = vec![Complex64::new(0.0, 0.0); N];
        for i in 0..N {
            let mut acc = Complex64::new(0.0, 0.0);
            for j in 0..N {
                acc += a[i + j * N] * x[j];
            }
            residual[i] = acc - b[i];
        }
        let nr: f64 = residual.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
        let nb: f64 = b.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
        let rel = nr / nb;
        assert!(rel < 1e-10, "relative residual {rel} >= 1e-10");
    }
}
