//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub use crate::blas::{Diag, Order, Side, Transpose, Uplo};
use crate::vector::{as_mut_ptr, as_ptr, check_equal_len, len, stride, Vector, VectorMut};

/// Modified matrix transformation (for the mathematical field `F`).
#[derive(Clone, Copy)]
pub enum H<F> {
    /// Specify that H is the matrix
    ///
    /// ⎧`h11`  `h12`⎫
    /// ⎩`h21`  `h22`⎭
    Full {
        h11: F,
        h21: F,
        h12: F,
        h22: F,
    },
    /// Specify that H is the matrix
    ///
    /// ⎧1.0  `h12`⎫
    /// ⎩`h21`  1.0⎭
    OffDiag {
        h21: F,
        h12: F,
    },
    /// Specify that H is the matrix
    ///
    /// ⎧`h11`   1.0⎫
    /// ⎩-1.0  `h22`⎭
    Diag {
        h11: F,
        h22: F,
    },
    Id,
}

/// `f32` vectors.
pub mod s {
    use super::*;

    // Level 1

    /// Return the sum of `alpha` and the dot product of `x` and `y`.
    #[doc(alias = "cblas_sdsdot")]
    pub fn sdot<T: Vector<f32> + ?Sized>(alpha: f32, x: &T, y: &T) -> f32 {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        unsafe { sys::cblas_sdsdot(len(x), alpha, as_ptr(x), stride(x), as_ptr(y), stride(y)) }
    }

    /// Return the dot product of `x` and `y`.
    #[doc(alias = "cblas_dsdot")]
    pub fn ddot<T: Vector<f32> + ?Sized>(x: &T, y: &T) -> f64 {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        unsafe { sys::cblas_dsdot(len(x), as_ptr(x), stride(x), as_ptr(y), stride(y)) }
    }
    /// Return the dot product of `x` and `y`.
    #[doc(alias = "cblas_sdot")]
    pub fn dot<T: Vector<f32> + ?Sized>(x: &T, y: &T) -> f32 {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        unsafe { sys::cblas_sdot(len(x), as_ptr(x), stride(x), as_ptr(y), stride(y)) }
    }

    /// Return the Euclidean norm of `x`.
    #[doc(alias = "cblas_snrm2")]
    pub fn nrm2<T: Vector<f32> + ?Sized>(x: &T) -> f32 {
        unsafe { sys::cblas_snrm2(len(x), as_ptr(x), stride(x)) }
    }

    /// Return the sum of the absolute values of the elements of `x`
    /// (i.e., its L¹-norm).
    #[doc(alias = "cblas_sasum")]
    pub fn asum<T: Vector<f32> + ?Sized>(x: &T) -> f32 {
        unsafe { sys::cblas_sasum(len(x), as_ptr(x), stride(x)) }
    }

    /// Return the index of the element with maximum absolute value.
    #[doc(alias = "cblas_isamax")]
    pub fn iamax<T: Vector<f32> + ?Sized>(x: &T) -> usize {
        unsafe { sys::cblas_isamax(len(x), as_ptr(x), stride(x)) }
    }

    /// Swap vectors `x` and `y`.
    #[doc(alias = "cblas_sswap")]
    pub fn swap<T1, T2>(x: &mut T1, y: &mut T2)
    where
        T1: VectorMut<f32> + ?Sized,
        T2: VectorMut<f32> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe { sys::cblas_sswap(len(x), as_mut_ptr(x), stride(x), as_mut_ptr(y), stride(y)) }
    }

    /// Copy the content of `x` into `y`.
    #[doc(alias = "cblas_scopy")]
    pub fn copy<T1, T2>(x: &T1, y: &mut T2)
    where
        T1: Vector<f32> + ?Sized,
        T2: VectorMut<f32> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe { sys::cblas_scopy(len(x), as_ptr(x), stride(x), as_mut_ptr(y), stride(y)) }
    }

    /// `y` := `alpha` * `x` + `y`.
    #[doc(alias = "cblas_saxpy")]
    pub fn axpy<T1, T2>(alpha: f32, x: &T1, y: &mut T2)
    where
        T1: Vector<f32> + ?Sized,
        T2: VectorMut<f32> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_saxpy(
                len(x),
                alpha,
                as_ptr(x),
                stride(x),
                as_mut_ptr(y),
                stride(y),
            )
        }
    }

    /// Given the Cartesian coordinates (`a`, `b`), returns
    /// (c, s, r, z) such that
    ///
    /// ⎧c  s⎫ ⎧a⎫ = ⎧r⎫
    /// ⎩s  c⎭ ⎩b⎭   ⎩0⎭
    ///
    /// The value of z is defined such that if |`a`| > |`b`|, z is s;
    /// otherwise if c ≠ 0, z is 1/c; otherwise z is 1.
    #[doc(alias = "cblas_srotg")]
    pub fn rotg(a: f32, b: f32) -> (f32, f32, f32, f32) {
        let mut c = f32::NAN;
        let mut s = f32::NAN;
        let mut r = a;
        let mut z = b;
        unsafe {
            sys::cblas_srotg(
                &mut r as *mut _,
                &mut z as *mut _,
                &mut c as *mut _,
                &mut s as *mut _,
            )
        }
        (c, s, r, z)
    }

    /// Given Cartesian coordinates (`x1`, `x2`), return the
    /// transformation matrix H that zeros the second component or the
    /// vector (`x1` √`d1`, `x2` √`d2`):
    ///
    /// H ⎧`x1` √`d1`⎫ = ⎧y1⎫
    ///   ⎩`x2` √`d2`⎭   ⎩0.⎭
    ///
    /// The second component of the return value is `y1`.
    #[doc(alias = "cblas_srotmg")]
    pub fn rotmg(mut d1: f32, mut d2: f32, mut x1: f32, x2: f32) -> (H<f32>, f32) {
        let mut h: [f32; 5] = [0.; 5];
        unsafe {
            sys::cblas_srotmg(
                &mut d1 as *mut _,
                &mut d2 as *mut _,
                &mut x1 as *mut _,
                x2,
                &mut h as *mut _,
            )
        }
        let h = if h[0] == -1.0 {
            H::Full {
                h11: h[1],
                h21: h[2],
                h12: h[3],
                h22: h[4],
            }
        } else if h[0] == 0.0 {
            H::OffDiag {
                h21: h[2],
                h12: h[3],
            }
        } else if h[0] == 1.0 {
            H::Diag {
                h11: h[1],
                h22: h[4],
            }
        } else if h[0] == -2.0 {
            H::Id
        } else {
            unreachable!("srotmg: incorrect flag value")
        };
        (h, x1)
    }

    /// Apply plane rotation.  More specifically, perform the
    /// following transformation in place :
    ///
    /// ⎧`x`ᵢ⎫ = ⎧`c`  `s`⎫ ⎧`x`ᵢ⎫
    /// ⎩`y`ᵢ⎭   ⎩-`s` `c`⎭ ⎩`y`ᵢ⎭
    ///
    /// for all indices i.
    #[doc(alias = "cblas_srot")]
    pub fn rot<T>(x: &mut T, y: &mut T, c: f32, s: f32)
    where
        T: VectorMut<f32> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_srot(
                len(x),
                as_mut_ptr(x),
                stride(x),
                as_mut_ptr(y),
                stride(y),
                c,
                s,
            )
        }
    }

    /// Apply the matrix rotation `h` to `x`, `y`.
    ///
    /// ⎧`x`ᵢ⎫ = `h` ⎧`x`ᵢ⎫
    /// ⎩`y`ᵢ⎭       ⎩`y`ᵢ⎭
    ///
    /// for all indices i.
    #[doc(alias = "cblas_srotm")]
    pub fn rotm<T>(x: &mut T, y: &mut T, h: H<f32>)
    where
        T: VectorMut<f32> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        let p = match h {
            H::Full { h11, h21, h12, h22 } => [-1.0, h11, h21, h12, h22],
            H::OffDiag { h21, h12 } => [0.0, 1., h21, h12, 1.],
            H::Diag { h11, h22 } => [1.0, h11, -1., 1., h22],
            H::Id => [-2.0, 1., 0., 0., 1.],
        };
        unsafe {
            sys::cblas_srotm(
                len(x),
                as_mut_ptr(x),
                stride(x),
                as_mut_ptr(y),
                stride(y),
                &p as *const _,
            )
        }
    }

    /// Multiply each element of `x` by `alpha`.
    #[doc(alias = "cblas_sscal")]
    pub fn scal<T>(alpha: f32, x: &mut T)
    where
        T: VectorMut<f32> + ?Sized,
    {
        unsafe { sys::cblas_sscal(len(x), alpha, as_mut_ptr(x), stride(x)) }
    }

    // Level 2

    /// Multiplies a matrix and a vector.
    ///
    /// * order : Whether matrices are row major order (C-Style) for column major order (Fortran-style). One of enum CblasRowMajor or CblasColMajor
    /// * transA :  Whether to transpose matrix A. One of enum CblasNoTrans, CBlasTrans.
    /// * M : Rows in matrix A
    /// * N : Columns in matrix A
    /// * alpha : scalar factor for (sigma * op(A) * x)
    /// * A : matrix A
    /// * lda : The size of the first dimension of matrix A
    /// * X : vector X
    /// * incx : use every incX'th element of X
    /// * beta : scalar factor y
    /// * Y : vector Y
    /// * incy : use every incY'th element of Y
    ///
    /// For parameter lda, if you are passing a matrix `A[m][n]`, the value of parameter lda should be m.
    #[doc(alias = "cblas_sgemv")]
    pub fn gemv(
        order: Order,
        transA: Transpose,
        M: i32,
        N: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        X: &[f32],
        incx: i32,
        beta: f32,
        Y: &mut [f32],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_sgemv(
                order.into(),
                transA.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                X.as_ptr(),
                incx,
                beta,
                Y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_sgbmv")]
    pub fn gbmv(
        order: Order,
        transA: Transpose,
        M: i32,
        N: i32,
        KL: i32,
        KU: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        X: &[f32],
        incx: i32,
        beta: f32,
        Y: &mut [f32],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_sgbmv(
                order.into(),
                transA.into(),
                M,
                N,
                KL,
                KU,
                alpha,
                A.as_ptr(),
                lda,
                X.as_ptr(),
                incx,
                beta,
                Y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_strmv")]
    pub fn trmv(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        A: &[f32],
        lda: i32,
        X: &mut [f32],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_strmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_stbmv")]
    pub fn tbmv(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        K: i32,
        A: &[f32],
        lda: i32,
        X: &mut [f32],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_stbmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_stpmv")]
    pub fn tpmv(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        Ap: &[f32],
        X: &mut [f32],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_stpmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr(),
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_strsv")]
    pub fn trsv(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        A: &[f32],
        lda: i32,
        X: &mut [f32],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_strsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_stbsv")]
    pub fn tbsv(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        K: i32,
        A: &[f32],
        lda: i32,
        X: &mut [f32],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_stbsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_stpsv")]
    pub fn tpsv(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        Ap: &[f32],
        X: &mut [f32],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_stpsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr(),
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ssymv")]
    pub fn symv(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        x: &[f32],
        incx: i32,
        beta: f32,
        y: &mut [f32],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_ssymv(
                order.into(),
                uplo.into(),
                N,
                alpha,
                A.as_ptr(),
                lda,
                x.as_ptr(),
                incx,
                beta,
                y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_ssbmv")]
    pub fn sbmv(
        order: Order,
        uplo: Uplo,
        N: i32,
        K: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        x: &[f32],
        incx: i32,
        beta: f32,
        y: &mut [f32],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_ssbmv(
                order.into(),
                uplo.into(),
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                x.as_ptr(),
                incx,
                beta,
                y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_sspmv")]
    pub fn spmv(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: f32,
        Ap: &[f32],
        x: &[f32],
        incx: i32,
        beta: f32,
        y: &mut [f32],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_sspmv(
                order.into(),
                uplo.into(),
                N,
                alpha,
                Ap.as_ptr(),
                x.as_ptr(),
                incx,
                beta,
                y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_sger")]
    pub fn ger(
        order: Order,
        M: i32,
        N: i32,
        alpha: f32,
        x: &[f32],
        incx: i32,
        y: &[f32],
        incy: i32,
        A: &mut [f32],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_sger(
                order.into(),
                M,
                N,
                alpha,
                x.as_ptr(),
                incx,
                y.as_ptr(),
                incy,
                A.as_mut_ptr(),
                lda,
            )
        }
    }

    #[doc(alias = "cblas_ssyr")]
    pub fn syr(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: f32,
        x: &[f32],
        incx: i32,
        A: &mut [f32],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_ssyr(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                A.as_mut_ptr(),
                lda,
            )
        }
    }

    #[doc(alias = "cblas_sspr")]
    pub fn spr(order: Order, uplo: Uplo, N: i32, alpha: f32, x: &[f32], incx: i32, Ap: &mut [f32]) {
        unsafe {
            sys::cblas_sspr(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                Ap.as_mut_ptr(),
            )
        }
    }

    #[doc(alias = "cblas_ssyr2")]
    pub fn syr2(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: f32,
        x: &[f32],
        incx: i32,
        y: &[f32],
        incy: i32,
        A: &mut [f32],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_ssyr2(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                y.as_ptr(),
                incy,
                A.as_mut_ptr(),
                lda,
            )
        }
    }

    #[doc(alias = "cblas_sspr2")]
    pub fn spr2(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: f32,
        x: &[f32],
        incx: i32,
        y: &[f32],
        incy: i32,
        A: &mut [f32],
    ) {
        unsafe {
            sys::cblas_sspr2(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                y.as_ptr(),
                incy,
                A.as_mut_ptr(),
            )
        }
    }

    // Level 3

    /// General crate::types::Matrix-MatrixF64 multiplication for single precision float.
    ///
    /// __Parameters:__
    ///
    /// * order : Whether matrices are row major order (C-Style) for column major order (Fortran-style). One of enum CblasRowMajor or CblasColMajor.
    /// * transA : Whether to transpose matrix A. One of enum CblasNoTrans, CBlasTrans, CBlasConjTrans.
    /// * transB : Whether to transpose matrix B. One of enum CblasNoTrans, CBlasTrans, CBlasConjTrans.
    /// * M : Rows in matrices A and C
    /// * N : Columns in Matrices B and C
    /// * K : Columns in matrix A and Rows in matrix B
    /// * alpha : scalar factor for op(A)op(B)
    /// * A : matrix A
    /// * lda : The size of the first dimension of matrix A
    /// * B : matrix B
    /// * ldb : The size of the first dimension of matrix B
    /// * beta : scalar factor for C
    /// * C : matrix C
    /// * ldc : The size of the first dimension of matrix C
    ///
    /// For parameters lda, ldb, and ldc, if you are passing a matrix `D[m][n]`, the value of parameter lda, ldb, or ldc should be m.
    #[doc(alias = "cblas_sgemm")]
    pub fn gemm(
        order: Order,
        transA: Transpose,
        transB: Transpose,
        M: i32,
        N: i32,
        K: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        B: &[f32],
        ldb: i32,
        beta: f32,
        C: &mut [f32],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_sgemm(
                order.into(),
                transA.into(),
                transB.into(),
                M,
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                B.as_ptr(),
                ldb,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    /// Symmetric crate::types::Matrix-MatrixF64 multiplication for single precision float.
    ///
    /// __Parameters:__
    ///
    /// * order : Whether matrices are row major order (C-Style) for column major order (Fortran-style). One of enum CblasRowMajor or CblasColMajor.
    /// * side : If CBlasSideLeft, perform (sigma(A)(B) + beta C). If CBlasSideRight, perform (sigma (B)(A) + beta C)
    /// * uplo : Indicates whether to use the upper (CBlasUpper) or lower (CBlasLower) triangle of matrix A
    /// * M : Rows in matrices A and C
    /// * N : Columns in Matrices B and C
    /// * alpha : scalar factor for op(A)op(B)
    /// * A : matrix A
    /// * lda : The size of the first dimension of matrix A
    /// * B : matrix B
    /// * ldb : The size of the first dimension of matrix B
    /// * beta : scalar factor for C
    /// * C : matrix C
    /// * ldc : The size of the first dimension of matrix C
    #[doc(alias = "cblas_ssymm")]
    pub fn symm(
        order: Order,
        side: Side,
        uplo: Uplo,
        M: i32,
        N: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        B: &[f32],
        ldb: i32,
        beta: f32,
        C: &mut [f32],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_ssymm(
                order.into(),
                side.into(),
                uplo.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                B.as_ptr(),
                ldb,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_ssyrk")]
    pub fn syrk(
        order: Order,
        uplo: Uplo,
        trans: Transpose,
        N: i32,
        K: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        beta: f32,
        C: &mut [f32],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_ssyrk(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_ssyr2k")]
    pub fn syr2k(
        order: Order,
        uplo: Uplo,
        trans: Transpose,
        N: i32,
        K: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        B: &[f32],
        ldb: i32,
        beta: f32,
        C: &mut [f32],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_ssyr2k(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                B.as_ptr(),
                ldb,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_strmm")]
    pub fn trmm(
        order: Order,
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        M: i32,
        N: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        B: &mut [f32],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_strmm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                B.as_mut_ptr(),
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_strsm")]
    pub fn trsm(
        order: Order,
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        M: i32,
        N: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        B: &mut [f32],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_strsm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                B.as_mut_ptr(),
                ldb,
            )
        }
    }
}

/// `f64` vectors.
pub mod d {
    use super::*;

    // Level 1

    /// Return the dot product of `x` and `y`.
    #[doc(alias = "cblas_ddot")]
    pub fn dot<T: Vector<f64> + ?Sized>(x: &T, y: &T) -> f64 {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        unsafe { sys::cblas_ddot(len(x), as_ptr(x), stride(x), as_ptr(y), stride(y)) }
    }
    /// Return the Euclidean norm of `x`.
    #[doc(alias = "cblas_dnrm2")]
    pub fn nrm2<T: Vector<f64> + ?Sized>(x: &T) -> f64 {
        unsafe { sys::cblas_dnrm2(len(x), as_ptr(x), stride(x)) }
    }

    /// Return the sum of the absolute values of the elements of `x`
    /// (i.e., its L¹-norm).
    #[doc(alias = "cblas_dasum")]
    pub fn asum<T: Vector<f64> + ?Sized>(x: &T) -> f64 {
        unsafe { sys::cblas_dasum(len(x), as_ptr(x), stride(x)) }
    }

    /// Return the index of the element with maximum absolute value.
    #[doc(alias = "cblas_idamax")]
    pub fn iamax<T: Vector<f64> + ?Sized>(x: &T) -> usize {
        unsafe { sys::cblas_idamax(len(x), as_ptr(x), stride(x)) }
    }

    /// Swap vectors `x` and `y`.
    #[doc(alias = "cblas_dswap")]
    pub fn swap<T1, T2>(x: &mut T1, y: &mut T2)
    where
        T1: VectorMut<f64> + ?Sized,
        T2: VectorMut<f64> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe { sys::cblas_dswap(len(x), as_mut_ptr(x), stride(x), as_mut_ptr(y), stride(y)) }
    }

    /// Copy the content of `x` into `y`.
    #[doc(alias = "cblas_dcopy")]
    pub fn copy<T1, T2>(x: &T1, y: &mut T2)
    where
        T1: Vector<f64> + ?Sized,
        T2: VectorMut<f64> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe { sys::cblas_dcopy(len(x), as_ptr(x), stride(x), as_mut_ptr(y), stride(y)) }
    }

    /// `y` := `alpha` * `x` + `y`.
    #[doc(alias = "cblas_daxpy")]
    pub fn axpy<T1, T2>(alpha: f64, x: &T1, y: &mut T2)
    where
        T1: Vector<f64> + ?Sized,
        T2: VectorMut<f64> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_daxpy(
                len(x),
                alpha,
                as_ptr(x),
                stride(x),
                as_mut_ptr(y),
                stride(y),
            )
        }
    }

    /// Given the Cartesian coordinates (`a`, `b`), returns
    /// (c, s, r, z) such that
    ///
    /// ⎧c  s⎫ ⎧a⎫ = ⎧r⎫
    /// ⎩s  c⎭ ⎩b⎭   ⎩0⎭
    ///
    /// The value of z is defined such that if |`a`| > |`b`|, z is s;
    /// otherwise if c ≠ 0, z is 1/c; otherwise z is 1.
    #[doc(alias = "cblas_drotg")]
    pub fn rotg(a: f64, b: f64) -> (f64, f64, f64, f64) {
        let mut c = f64::NAN;
        let mut s = f64::NAN;
        let mut r = a;
        let mut z = b;
        unsafe {
            sys::cblas_drotg(
                &mut r as *mut _,
                &mut z as *mut _,
                &mut c as *mut _,
                &mut s as *mut _,
            )
        }
        (c, s, r, z)
    }

    /// Given Cartesian coordinates (`x1`, `x2`), return the
    /// transformation matrix H that zeros the second component or the
    /// vector (`x1` √`d1`, `x2` √`d2`):
    ///
    /// H ⎧`x1` √`d1`⎫ = ⎧y1⎫
    ///   ⎩`x2` √`d2`⎭   ⎩0.⎭
    ///
    /// The second component of the return value is `y1`.
    #[doc(alias = "cblas_drotmg")]
    pub fn rotmg(mut d1: f64, mut d2: f64, mut x1: f64, x2: f64) -> (H<f64>, f64) {
        let mut h: [f64; 5] = [0.; 5];
        unsafe {
            sys::cblas_drotmg(
                &mut d1 as *mut _,
                &mut d2 as *mut _,
                &mut x1 as *mut _,
                x2,
                &mut h as *mut _,
            )
        }
        let h = if h[0] == -1.0 {
            H::Full {
                h11: h[1],
                h21: h[2],
                h12: h[3],
                h22: h[4],
            }
        } else if h[0] == 0.0 {
            H::OffDiag {
                h21: h[2],
                h12: h[3],
            }
        } else if h[0] == 1.0 {
            H::Diag {
                h11: h[1],
                h22: h[4],
            }
        } else if h[0] == -2.0 {
            H::Id
        } else {
            unreachable!("srotmg: incorrect flag value")
        };
        (h, x1)
    }

    /// Apply plane rotation.  More specifically, perform the
    /// following transformation in place :
    ///
    /// ⎧`x`ᵢ⎫ = ⎧`c`  `s`⎫ ⎧`x`ᵢ⎫
    /// ⎩`y`ᵢ⎭   ⎩-`s` `c`⎭ ⎩`y`ᵢ⎭
    ///
    /// for all indices i.
    #[doc(alias = "cblas_drot")]
    pub fn rot<T>(x: &mut T, y: &mut T, c: f64, s: f64)
    where
        T: VectorMut<f64> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_drot(
                len(x),
                as_mut_ptr(x),
                stride(x),
                as_mut_ptr(y),
                stride(y),
                c,
                s,
            )
        }
    }

    /// Apply the matrix rotation `h` to `x`, `y`.
    ///
    /// ⎧`x`ᵢ⎫ = `h` ⎧`x`ᵢ⎫
    /// ⎩`y`ᵢ⎭       ⎩`y`ᵢ⎭
    ///
    /// for all indices i.
    #[doc(alias = "cblas_drotm")]
    pub fn rotm<T>(x: &mut T, y: &mut T, h: H<f64>)
    where
        T: VectorMut<f64> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        let p = match h {
            H::Full { h11, h21, h12, h22 } => [-1.0, h11, h21, h12, h22],
            H::OffDiag { h21, h12 } => [0.0, 1., h21, h12, 1.],
            H::Diag { h11, h22 } => [1.0, h11, -1., 1., h22],
            H::Id => [-2.0, 1., 0., 0., 1.],
        };
        unsafe {
            sys::cblas_drotm(
                len(x),
                as_mut_ptr(x),
                stride(x),
                as_mut_ptr(y),
                stride(y),
                &p as *const _,
            )
        }
    }

    /// Multiply each element of `x` by `alpha`.
    #[doc(alias = "cblas_dscal")]
    pub fn scal<T>(alpha: f64, x: &mut T)
    where
        T: VectorMut<f64> + ?Sized,
    {
        unsafe { sys::cblas_dscal(len(x), alpha, as_mut_ptr(x), stride(x)) }
    }

    // Level 2

    #[doc(alias = "cblas_dgemv")]
    pub fn gemv(
        order: Order,
        transA: Transpose,
        M: i32,
        N: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        X: &[f64],
        incx: i32,
        beta: f64,
        Y: &mut [f64],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_dgemv(
                order.into(),
                transA.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                X.as_ptr(),
                incx,
                beta,
                Y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_dgbmv")]
    pub fn gbmv(
        order: Order,
        transA: Transpose,
        M: i32,
        N: i32,
        KL: i32,
        KU: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        X: &[f64],
        incx: i32,
        beta: f64,
        Y: &mut [f64],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_dgbmv(
                order.into(),
                transA.into(),
                M,
                N,
                KL,
                KU,
                alpha,
                A.as_ptr(),
                lda,
                X.as_ptr(),
                incx,
                beta,
                Y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_dtrmv")]
    pub fn trmv(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        A: &[f64],
        lda: i32,
        X: &mut [f64],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_dtrmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_dtbmv")]
    pub fn tbmv(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        K: i32,
        A: &[f64],
        lda: i32,
        X: &mut [f64],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_dtbmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_dtpmv")]
    pub fn tpmv(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        Ap: &[f64],
        X: &mut [f64],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_dtpmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr(),
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_dtrsv")]
    pub fn trsv(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        A: &[f64],
        lda: i32,
        X: &mut [f64],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_dtrsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_dtbsv")]
    pub fn tbsv(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        K: i32,
        A: &[f64],
        lda: i32,
        X: &mut [f64],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_dtbsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_dtpsv")]
    pub fn tpsv(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        Ap: &[f64],
        X: &mut [f64],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_dtpsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr(),
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_dsymv")]
    pub fn symv(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        x: &[f64],
        incx: i32,
        beta: f64,
        y: &mut [f64],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_dsymv(
                order.into(),
                uplo.into(),
                N,
                alpha,
                A.as_ptr(),
                lda,
                x.as_ptr(),
                incx,
                beta,
                y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_dsbmv")]
    pub fn sbmv(
        order: Order,
        uplo: Uplo,
        N: i32,
        K: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        x: &[f64],
        incx: i32,
        beta: f64,
        y: &mut [f64],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_dsbmv(
                order.into(),
                uplo.into(),
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                x.as_ptr(),
                incx,
                beta,
                y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_dspmv")]
    pub fn spmv(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: f64,
        Ap: &[f64],
        x: &[f64],
        incx: i32,
        beta: f64,
        y: &mut [f64],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_dspmv(
                order.into(),
                uplo.into(),
                N,
                alpha,
                Ap.as_ptr(),
                x.as_ptr(),
                incx,
                beta,
                y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_dger")]
    pub fn ger(
        order: Order,
        M: i32,
        N: i32,
        alpha: f64,
        x: &[f64],
        incx: i32,
        y: &[f64],
        incy: i32,
        A: &mut [f64],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_dger(
                order.into(),
                M,
                N,
                alpha,
                x.as_ptr(),
                incx,
                y.as_ptr(),
                incy,
                A.as_mut_ptr(),
                lda,
            )
        }
    }

    #[doc(alias = "cblas_dsyr")]
    pub fn syr(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: f64,
        x: &[f64],
        incx: i32,
        A: &mut [f64],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_dsyr(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                A.as_mut_ptr(),
                lda,
            )
        }
    }

    #[doc(alias = "cblas_dspr")]
    pub fn spr(order: Order, uplo: Uplo, N: i32, alpha: f64, x: &[f64], incx: i32, Ap: &mut [f64]) {
        unsafe {
            sys::cblas_dspr(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                Ap.as_mut_ptr(),
            )
        }
    }

    #[doc(alias = "cblas_dsyr2")]
    pub fn syr2(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: f64,
        x: &[f64],
        incx: i32,
        y: &[f64],
        incy: i32,
        A: &mut [f64],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_dsyr2(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                y.as_ptr(),
                incy,
                A.as_mut_ptr(),
                lda,
            )
        }
    }

    #[doc(alias = "cblas_dspr2")]
    pub fn spr2(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: f64,
        x: &[f64],
        incx: i32,
        y: &[f64],
        incy: i32,
        A: &mut [f64],
    ) {
        unsafe {
            sys::cblas_dspr2(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                y.as_ptr(),
                incy,
                A.as_mut_ptr(),
            )
        }
    }

    // Level 3

    #[doc(alias = "cblas_dgemm")]
    pub fn gemm(
        order: Order,
        transA: Transpose,
        transB: Transpose,
        M: i32,
        N: i32,
        K: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        B: &[f64],
        ldb: i32,
        beta: f64,
        C: &mut [f64],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_dgemm(
                order.into(),
                transA.into(),
                transB.into(),
                M,
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                B.as_ptr(),
                ldb,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_dsymm")]
    pub fn symm(
        order: Order,
        side: Side,
        uplo: Uplo,
        M: i32,
        N: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        B: &[f64],
        ldb: i32,
        beta: f64,
        C: &mut [f64],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_dsymm(
                order.into(),
                side.into(),
                uplo.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                B.as_ptr(),
                ldb,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_dsyrk")]
    pub fn syrk(
        order: Order,
        uplo: Uplo,
        trans: Transpose,
        N: i32,
        K: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        beta: f64,
        C: &mut [f64],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_dsyrk(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_dsyr2k")]
    pub fn syr2k(
        order: Order,
        uplo: Uplo,
        trans: Transpose,
        N: i32,
        K: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        B: &[f64],
        ldb: i32,
        beta: f64,
        C: &mut [f64],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_dsyr2k(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                B.as_ptr(),
                ldb,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_dtrmm")]
    pub fn trmm(
        order: Order,
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        M: i32,
        N: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        B: &mut [f64],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_dtrmm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                B.as_mut_ptr(),
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_dtrsm")]
    pub fn trsm(
        order: Order,
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        M: i32,
        N: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        B: &mut [f64],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_dtrsm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                B.as_mut_ptr(),
                ldb,
            )
        }
    }
}

/// `Complex<f32>` vectors.
#[cfg(feature = "complex")]
#[cfg_attr(docsrs, doc(cfg(feature = "complex")))]
pub mod c {
    use super::*;
    use num_complex::Complex;

    // Level 1

    /// Return the unconjugated dot product between `x` and `y`, that
    /// is $∑ x_i y_i$.
    ///
    /// # Example
    ///
    /// ```
    /// use num_complex::Complex;
    /// use rgsl::cblas::c;
    /// let x = [Complex::new(1., 1.), Complex::new(2., 1.)];
    /// assert_eq!(c::dotu(&x, &x), Complex::new(3., 6.))
    /// ```
    #[doc(alias = "cblas_cdotu_sub")]
    pub fn dotu<T>(x: &T, y: &T) -> Complex<f32>
    where
        T: Vector<Complex<f32>> + ?Sized,
    {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        let mut dotu: Complex<f32> = Complex::new(0., 0.);
        unsafe {
            sys::cblas_cdotu_sub(
                len(x),
                as_ptr(x) as *const _,
                stride(x),
                as_ptr(y) as *const _,
                stride(y),
                &mut dotu as *mut Complex<f32> as *mut _,
            )
        }
        dotu
    }

    /// Return the (conjugated) dot product between `x` and `y`, that
    /// is ∑ x̅ᵢ yᵢ.
    ///
    /// # Example
    ///
    /// ```
    /// use num_complex::Complex;
    /// use rgsl::cblas::c;
    /// let x = [Complex::new(1., 1.), Complex::new(2., 1.)];
    /// assert_eq!(c::dot(&x, &x), Complex::new(7., 0.))
    /// ```
    #[doc(alias = "cblas_cdotc_sub")]
    pub fn dot<T>(x: &T, y: &T) -> Complex<f32>
    where
        T: Vector<Complex<f32>> + ?Sized,
    {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        let mut dotc: Complex<f32> = Complex::new(0., 0.);
        unsafe {
            sys::cblas_cdotc_sub(
                len(x),
                as_ptr(x) as *const _,
                stride(x),
                as_ptr(y) as *const _,
                stride(y),
                &mut dotc as *mut Complex<f32> as *mut _,
            )
        }
        dotc
    }

    /// Return the Euclidean norm of `x`.
    ///
    /// # Example
    ///
    /// ```
    /// use num_complex::Complex;
    /// use rgsl::cblas::c;
    /// let x = [Complex::new(1., 1.), Complex::new(2., 1.)];
    /// assert_eq!(c::nrm2(&x), 7f32.sqrt())
    /// ```
    #[doc(alias = "cblas_scnrm2")]
    pub fn nrm2<T: Vector<Complex<f32>> + ?Sized>(x: &T) -> f32 {
        unsafe { sys::cblas_scnrm2(len(x), as_ptr(x) as *const _, stride(x)) }
    }

    /// Return the sum of the modulus of the elements of `x`
    /// (i.e., its L¹-norm).
    #[doc(alias = "cblas_scasum")]
    pub fn asum<T: Vector<Complex<f32>> + ?Sized>(x: &T) -> f32 {
        unsafe { sys::cblas_scasum(len(x), as_ptr(x) as *const _, stride(x)) }
    }

    /// Return the index of the element with maximum modulus.
    #[doc(alias = "cblas_icamax")]
    pub fn iamax<T: Vector<Complex<f32>> + ?Sized>(x: &T) -> usize {
        unsafe { sys::cblas_icamax(len(x), as_ptr(x) as *const _, stride(x)) }
    }

    /// Swap vectors `x` and `y`.
    #[doc(alias = "cblas_cswap")]
    pub fn swap<T1, T2>(x: &mut T1, y: &mut T2)
    where
        T1: VectorMut<Complex<f32>> + ?Sized,
        T2: VectorMut<Complex<f32>> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_cswap(
                len(x),
                as_mut_ptr(x) as *mut _,
                stride(x),
                as_mut_ptr(y) as *mut _,
                stride(y),
            )
        }
    }

    /// Copy the content of `x` into `y`.
    #[doc(alias = "cblas_ccopy")]
    pub fn copy<T1, T2>(x: &T1, y: &mut T2)
    where
        T1: Vector<Complex<f32>> + ?Sized,
        T2: VectorMut<Complex<f32>> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_ccopy(
                len(x),
                as_ptr(x) as *const _,
                stride(x),
                as_mut_ptr(y) as *mut _,
                stride(y),
            )
        }
    }

    /// `y` := `alpha` * `x` + `y`.
    #[doc(alias = "cblas_caxpy")]
    pub fn axpy<T1, T2>(alpha: &Complex<f32>, x: &T1, y: &mut T2)
    where
        T1: Vector<Complex<f32>> + ?Sized,
        T2: VectorMut<Complex<f32>> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_caxpy(
                len(x),
                alpha as *const Complex<f32> as *const _,
                as_ptr(x) as *const _,
                stride(x),
                as_mut_ptr(y) as *mut _,
                stride(y),
            )
        }
    }

    /// Multiply each element of `x` by `alpha`.
    #[doc(alias = "cblas_cscal")]
    pub fn scal<T>(alpha: &Complex<f32>, x: &mut T)
    where
        T: VectorMut<Complex<f32>> + ?Sized,
    {
        unsafe {
            sys::cblas_cscal(
                len(x),
                alpha as *const Complex<f32> as *const _,
                as_mut_ptr(x) as *mut _,
                stride(x),
            )
        }
    }

    /// Multiply each element of `x` by `alpha`.
    #[doc(alias = "cblas_csscal")]
    pub fn rscal<T>(alpha: f32, x: &mut T)
    where
        T: VectorMut<Complex<f32>> + ?Sized,
    {
        unsafe { sys::cblas_csscal(len(x), alpha, as_mut_ptr(x) as *mut _, stride(x)) }
    }

    #[doc(alias = "cblas_cgemv")]
    pub fn gemv<T>(
        order: Order,
        transA: Transpose,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        X: &[T],
        incx: i32,
        beta: &[T],
        Y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_cgemv(
                order.into(),
                transA.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                X.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                Y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_cgbmv")]
    pub fn gbmv<T>(
        order: Order,
        transA: Transpose,
        M: i32,
        N: i32,
        KL: i32,
        KU: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        X: &[T],
        incx: i32,
        beta: &[T],
        Y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_cgbmv(
                order.into(),
                transA.into(),
                M,
                N,
                KL,
                KU,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                X.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                Y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_ctrmv")]
    pub fn trmv<T>(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ctrmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ctbmv")]
    pub fn tbmv<T>(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        K: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ctbmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ctpmv")]
    pub fn tpmv<T>(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        Ap: &[T],
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ctpmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr() as *const _,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ctrsv")]
    pub fn trsv<T>(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ctrsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ctbsv")]
    pub fn tbsv<T>(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        K: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ctbsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ctpsv")]
    pub fn tpsv<T>(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        Ap: &[T],
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ctpsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr() as *const _,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_chemv")]
    pub fn hemv<T>(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        x: &[T],
        incx: i32,
        beta: &[T],
        y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_chemv(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                x.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_chbmv")]
    pub fn hbmv<T>(
        order: Order,
        uplo: Uplo,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        x: &[T],
        incx: i32,
        beta: &[T],
        y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_chbmv(
                order.into(),
                uplo.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                x.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_chpmv")]
    pub fn hpmv<T>(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: &[T],
        Ap: &[T],
        x: &[T],
        incx: i32,
        beta: &[T],
        y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_chpmv(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                Ap.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_cgeru")]
    pub fn geru<T>(
        order: Order,
        M: i32,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[T],
        incy: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_cgeru(
                order.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_cgerc")]
    pub fn gerc<T>(
        order: Order,
        M: i32,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[T],
        incy: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_cgerc(
                order.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_cher")]
    pub fn her<T>(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: f32,
        x: &[T],
        incx: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_cher(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr() as *const _,
                incx,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_chpr")]
    pub fn hpr<T>(order: Order, uplo: Uplo, N: i32, alpha: f32, x: &[T], incx: i32, Ap: &mut [T]) {
        unsafe {
            sys::cblas_chpr(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr() as *const _,
                incx,
                Ap.as_mut_ptr() as *mut _,
            )
        }
    }

    #[doc(alias = "cblas_cher2")]
    pub fn her2<T>(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[T],
        incy: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_cher2(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_chpr2")]
    pub fn hpr2<T>(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[f64],
        incy: i32,
        Ap: &mut [f64],
    ) {
        unsafe {
            sys::cblas_chpr2(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                Ap.as_mut_ptr() as *mut _,
            )
        }
    }

    // Level 3

    #[doc(alias = "cblas_cgemm")]
    pub fn gemm<T>(
        order: Order,
        transA: Transpose,
        transB: Transpose,
        M: i32,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_cgemm(
                order.into(),
                transA.into(),
                transB.into(),
                M,
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_csymm")]
    pub fn symm<T>(
        order: Order,
        side: Side,
        uplo: Uplo,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_csymm(
                order.into(),
                side.into(),
                uplo.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_csyrk")]
    pub fn syrk<T>(
        order: Order,
        uplo: Uplo,
        trans: Transpose,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_csyrk(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_csyr2k")]
    pub fn syr2k<T>(
        order: Order,
        uplo: Uplo,
        trans: Transpose,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_csyr2k(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_ctrmm")]
    pub fn trmm<T>(
        order: Order,
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &mut [T],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_ctrmm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_mut_ptr() as *mut _,
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_ctrsm")]
    pub fn trsm<T>(
        order: Order,
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &mut [T],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_ctrsm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_mut_ptr() as *mut _,
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_chemm")]
    pub fn hemm<T>(
        order: Order,
        side: Side,
        uplo: Uplo,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_chemm(
                order.into(),
                side.into(),
                uplo.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_cherk")]
    pub fn herk<T>(
        order: Order,
        uplo: Uplo,
        trans: Transpose,
        N: i32,
        K: i32,
        alpha: f32,
        A: &[T],
        lda: i32,
        beta: f32,
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_cherk(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha,
                A.as_ptr() as *const _,
                lda,
                beta,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_cher2k")]
    pub fn her2k<T>(
        order: Order,
        uplo: Uplo,
        trans: Transpose,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: f32,
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_cher2k(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }
}

/// `Complex<f64>` vectors.
#[cfg(feature = "complex")]
#[cfg_attr(docsrs, doc(cfg(feature = "complex")))]
pub mod z {
    use super::*;
    use num_complex::Complex;

    // Level 1

    /// Return the unconjugated dot product between `x` and `y`, that
    /// is ∑ xᵢ yᵢ.
    ///
    /// # Example
    ///
    /// ```
    /// use num_complex::Complex;
    /// use rgsl::cblas::z;
    /// let x = [Complex::new(1., 1.), Complex::new(2., 1.)];
    /// assert_eq!(z::dotu(&x, &x), Complex::new(3., 6.))
    /// ```
    #[doc(alias = "cblas_zdotu_sub")]
    pub fn dotu<T>(x: &T, y: &T) -> Complex<f64>
    where
        T: Vector<Complex<f64>> + ?Sized,
    {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        let mut dotu: Complex<f64> = Complex::new(0., 0.);
        unsafe {
            sys::cblas_zdotu_sub(
                len(x),
                as_ptr(x) as *const _,
                stride(x),
                as_ptr(y) as *const _,
                stride(y),
                &mut dotu as *mut Complex<f64> as *mut _,
            )
        }
        dotu
    }

    /// Return the (conjugated) dot product between `x` and `y`, that
    /// is ∑ x̅ᵢ yᵢ.
    ///
    /// # Example
    ///
    /// ```
    /// use num_complex::Complex;
    /// use rgsl::cblas::z;
    /// let x = [Complex::new(1., 1.), Complex::new(2., 1.)];
    /// assert_eq!(z::dot(&x, &x), Complex::new(7., 0.))
    /// ```
    #[doc(alias = "cblas_zdotc_sub")]
    pub fn dot<T>(x: &T, y: &T) -> Complex<f64>
    where
        T: Vector<Complex<f64>> + ?Sized,
    {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        let mut dotc: Complex<f64> = Complex::new(0., 0.);
        unsafe {
            sys::cblas_zdotc_sub(
                len(x),
                as_ptr(x) as *const _,
                stride(x),
                as_ptr(y) as *const _,
                stride(y),
                &mut dotc as *mut Complex<f64> as *mut _,
            )
        }
        dotc
    }

    /// Return the Euclidean norm of `x`.
    ///
    /// # Example
    ///
    /// ```
    /// use num_complex::Complex;
    /// use rgsl::cblas::z;
    /// let x = [Complex::new(1., 1.), Complex::new(2., 1.)];
    /// assert_eq!(z::nrm2(&x), 7f64.sqrt())
    /// ```
    #[doc(alias = "cblas_dznrm2")]
    pub fn nrm2<T: Vector<Complex<f64>> + ?Sized>(x: &T) -> f64 {
        unsafe { sys::cblas_dznrm2(len(x), as_ptr(x) as *const _, stride(x)) }
    }

    /// Return the sum of the modulus of the elements of `x`
    /// (i.e., its L¹-norm).
    #[doc(alias = "cblas_dzasum")]
    pub fn asum<T: Vector<Complex<f64>> + ?Sized>(x: &T) -> f64 {
        unsafe { sys::cblas_dzasum(len(x), as_ptr(x) as *const _, stride(x)) }
    }

    /// Return the index of the element with maximum modulus.
    #[doc(alias = "cblas_izamax")]
    pub fn iamax<T: Vector<Complex<f64>> + ?Sized>(x: &T) -> usize {
        unsafe { sys::cblas_izamax(len(x), as_ptr(x) as *const _, stride(x)) }
    }

    /// Swap vectors `x` and `y`.
    #[doc(alias = "cblas_zswap")]
    pub fn swap<T1, T2>(x: &mut T1, y: &mut T2)
    where
        T1: VectorMut<Complex<f64>> + ?Sized,
        T2: VectorMut<Complex<f64>> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_zswap(
                len(x),
                as_mut_ptr(x) as *mut _,
                stride(x),
                as_mut_ptr(y) as *mut _,
                stride(y),
            )
        }
    }

    /// Copy the content of `x` into `y`.
    #[doc(alias = "cblas_zcopy")]
    pub fn copy<T1, T2>(x: &T1, y: &mut T2)
    where
        T1: Vector<Complex<f64>> + ?Sized,
        T2: VectorMut<Complex<f64>> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_zcopy(
                len(x),
                as_ptr(x) as *const _,
                stride(x),
                as_mut_ptr(y) as *mut _,
                stride(y),
            )
        }
    }

    /// `y` := `alpha` * `x` + `y`.
    #[doc(alias = "cblas_zaxpy")]
    pub fn axpy<T1, T2>(alpha: &Complex<f64>, x: &T1, y: &mut T2)
    where
        T1: Vector<Complex<f64>> + ?Sized,
        T2: VectorMut<Complex<f64>> + ?Sized,
    {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_zaxpy(
                len(x),
                alpha as *const Complex<f64> as *const _,
                as_ptr(x) as *const _,
                stride(x),
                as_mut_ptr(y) as *mut _,
                stride(y),
            )
        }
    }

    /// Multiply each element of `x` by `alpha`.
    #[doc(alias = "cblas_zscal")]
    pub fn scal<T>(alpha: &Complex<f64>, x: &mut T)
    where
        T: VectorMut<Complex<f64>> + ?Sized,
    {
        unsafe {
            sys::cblas_zscal(
                len(x),
                alpha as *const Complex<f64> as *const _,
                as_mut_ptr(x) as *mut _,
                stride(x),
            )
        }
    }

    #[cfg(feature = "complex")]
    /// Multiple each element of a matrix/vector by a constant.
    #[doc(alias = "cblas_zdscal")]
    pub fn rscal<T>(alpha: f64, x: &mut T)
    where
        T: VectorMut<Complex<f64>> + ?Sized,
    {
        unsafe { sys::cblas_zdscal(len(x), alpha, as_mut_ptr(x) as *mut _, stride(x)) }
    }

    // Level 2

    #[doc(alias = "cblas_zgemv")]
    pub fn gemv<T>(
        order: Order,
        transA: Transpose,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        X: &[T],
        incx: i32,
        beta: &[T],
        Y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_zgemv(
                order.into(),
                transA.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                X.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                Y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_zgbmv")]
    pub fn gbmv<T>(
        order: Order,
        transA: Transpose,
        M: i32,
        N: i32,
        KL: i32,
        KU: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        X: &[T],
        incx: i32,
        beta: &[T],
        Y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_zgbmv(
                order.into(),
                transA.into(),
                M,
                N,
                KL,
                KU,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                X.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                Y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_ztrmv")]
    pub fn trmv<T>(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ztrmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ztbmv")]
    pub fn tbmv<T>(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        K: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ztbmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ztpmv")]
    pub fn tpmv<T>(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        Ap: &[T],
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ztpmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr() as *const _,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ztrsv")]
    pub fn trsv<T>(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ztrsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ztbsv")]
    pub fn tbsv<T>(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        K: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ztbsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ztpsv")]
    pub fn tpsv<T>(
        order: Order,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        N: i32,
        Ap: &[T],
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ztpsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr() as *const _,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_zhemv")]
    pub fn hemv<T>(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        x: &[T],
        incx: i32,
        beta: &[T],
        y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_zhemv(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                x.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_zhbmv")]
    pub fn hbmv<T>(
        order: Order,
        uplo: Uplo,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        x: &[T],
        incx: i32,
        beta: &[T],
        y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_zhbmv(
                order.into(),
                uplo.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                x.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_zhpmv")]
    pub fn hpmv<T>(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: &[T],
        Ap: &[T],
        x: &[T],
        incx: i32,
        beta: &[T],
        y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_zhpmv(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                Ap.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_zgeru")]
    pub fn geru<T>(
        order: Order,
        M: i32,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[T],
        incy: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_zgeru(
                order.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_zgerc")]
    pub fn gerc<T>(
        order: Order,
        M: i32,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[T],
        incy: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_zgerc(
                order.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_zher")]
    pub fn her<T>(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: f64,
        x: &[T],
        incx: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_zher(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr() as *const _,
                incx,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_zhpr")]
    pub fn hpr<T>(order: Order, uplo: Uplo, N: i32, alpha: f64, x: &[T], incx: i32, Ap: &mut [T]) {
        unsafe {
            sys::cblas_zhpr(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr() as *const _,
                incx,
                Ap.as_mut_ptr() as *mut _,
            )
        }
    }

    #[doc(alias = "cblas_zher2")]
    pub fn her2<T>(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[T],
        incy: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_zher2(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_zhpr2")]
    pub fn hpr2<T>(
        order: Order,
        uplo: Uplo,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[f64],
        incy: i32,
        Ap: &mut [f64],
    ) {
        unsafe {
            sys::cblas_zhpr2(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                Ap.as_mut_ptr() as *mut _,
            )
        }
    }

    // Level 3

    #[doc(alias = "cblas_zgemm")]
    pub fn gemm<T>(
        order: Order,
        transA: Transpose,
        transB: Transpose,
        M: i32,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zgemm(
                order.into(),
                transA.into(),
                transB.into(),
                M,
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_zsymm")]
    pub fn symm<T>(
        order: Order,
        side: Side,
        uplo: Uplo,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zsymm(
                order.into(),
                side.into(),
                uplo.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_zsyrk")]
    pub fn syrk<T>(
        order: Order,
        uplo: Uplo,
        trans: Transpose,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zsyrk(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_zsyr2k")]
    pub fn syr2k<T>(
        order: Order,
        uplo: Uplo,
        trans: Transpose,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zsyr2k(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_ztrmm")]
    pub fn trmm<T>(
        order: Order,
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &mut [T],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_ztrmm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_mut_ptr() as *mut _,
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_ztrsm")]
    pub fn trsm<T>(
        order: Order,
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &mut [T],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_ztrsm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_mut_ptr() as *mut _,
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_zhemm")]
    pub fn hemm<T>(
        order: Order,
        side: Side,
        uplo: Uplo,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zhemm(
                order.into(),
                side.into(),
                uplo.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_zherk")]
    pub fn herk<T>(
        order: Order,
        uplo: Uplo,
        trans: Transpose,
        N: i32,
        K: i32,
        alpha: f64,
        A: &[T],
        lda: i32,
        beta: f64,
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zherk(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha,
                A.as_ptr() as *const _,
                lda,
                beta,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_zher2k")]
    pub fn her2k<T>(
        order: Order,
        uplo: Uplo,
        trans: Transpose,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: f64,
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zher2k(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }
}
