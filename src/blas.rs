//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum Transpose {
    NoTranspose,
    Transpose,
    ConjugateTranspose,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::CBLAS_TRANSPOSE> for Transpose {
    fn into(self) -> sys::CBLAS_TRANSPOSE {
        match self {
            Self::NoTranspose => sys::CBLAS_TRANSPOSE_CblasNoTrans,
            Self::Transpose => sys::CBLAS_TRANSPOSE_CblasTrans,
            Self::ConjugateTranspose => sys::CBLAS_TRANSPOSE_CblasConjTrans,
        }
    }
}

#[doc(hidden)]
impl From<sys::CBLAS_TRANSPOSE> for Transpose {
    fn from(v: sys::CBLAS_TRANSPOSE) -> Transpose {
        match v {
            sys::CBLAS_TRANSPOSE_CblasNoTrans => Self::NoTranspose,
            sys::CBLAS_TRANSPOSE_CblasTrans => Self::Transpose,
            sys::CBLAS_TRANSPOSE_CblasConjTrans => Self::ConjugateTranspose,
            _ => panic!("Unknown CblasTranspose value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum Uplo {
    Upper,
    Lower,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::CBLAS_UPLO> for Uplo {
    fn into(self) -> sys::CBLAS_UPLO {
        match self {
            Self::Upper => sys::CBLAS_UPLO_CblasUpper,
            Self::Lower => sys::CBLAS_UPLO_CblasLower,
        }
    }
}

#[doc(hidden)]
impl From<sys::CBLAS_UPLO> for Uplo {
    fn from(v: sys::CBLAS_UPLO) -> Uplo {
        match v {
            sys::CBLAS_UPLO_CblasUpper => Self::Upper,
            sys::CBLAS_UPLO_CblasLower => Self::Lower,
            _ => panic!("Unknown CblasUplo value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum Order {
    RowMajor,
    ColumnMajor,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::CBLAS_ORDER> for Order {
    fn into(self) -> sys::CBLAS_ORDER {
        match self {
            Self::RowMajor => sys::CBLAS_ORDER_CblasRowMajor,
            Self::ColumnMajor => sys::CBLAS_ORDER_CblasColMajor,
        }
    }
}

#[doc(hidden)]
impl From<sys::CBLAS_ORDER> for Order {
    fn from(v: sys::CBLAS_ORDER) -> Order {
        match v {
            sys::CBLAS_ORDER_CblasRowMajor => Self::RowMajor,
            sys::CBLAS_ORDER_CblasColMajor => Self::ColumnMajor,
            _ => panic!("Unknown CblasOrder value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum Side {
    Left,
    Right,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::CBLAS_SIDE> for Side {
    fn into(self) -> sys::CBLAS_SIDE {
        match self {
            Self::Left => sys::CBLAS_SIDE_CblasLeft,
            Self::Right => sys::CBLAS_SIDE_CblasRight,
        }
    }
}

#[doc(hidden)]
impl From<sys::CBLAS_SIDE> for Side {
    fn from(v: sys::CBLAS_SIDE) -> Side {
        match v {
            sys::CBLAS_SIDE_CblasLeft => Self::Left,
            sys::CBLAS_SIDE_CblasRight => Self::Right,
            _ => panic!("Unknown CblasSide value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum Diag {
    NonUnit,
    Unit,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::CBLAS_SIDE> for Diag {
    fn into(self) -> sys::CBLAS_SIDE {
        match self {
            Self::NonUnit => sys::CBLAS_DIAG_CblasNonUnit,
            Self::Unit => sys::CBLAS_DIAG_CblasUnit,
        }
    }
}

#[doc(hidden)]
impl From<sys::CBLAS_SIDE> for Diag {
    fn from(v: sys::CBLAS_SIDE) -> Diag {
        match v {
            sys::CBLAS_DIAG_CblasNonUnit => Self::NonUnit,
            sys::CBLAS_DIAG_CblasUnit => Self::Unit,
            _ => panic!("Unknown CblasDiag value"),
        }
    }
}

/// `f32` vectors.
pub mod s {
    use super::*;
    use crate::{ffi::FFI, Error, MatrixF32, VectorF32};

    // Level 1

    /// Return the sum `alpha` + `x`ᵀ `y` of the vectors `x` and
    /// `y`.
    #[doc(alias = "gsl_blas_sdsdot")]
    pub fn sdot(alpha: f32, x: &VectorF32, y: &VectorF32) -> Result<f32, Error> {
        let mut result = 0.;
        let ret = unsafe {
            sys::gsl_blas_sdsdot(alpha, x.unwrap_shared(), y.unwrap_shared(), &mut result)
        };
        Error::handle(ret, result)
    }

    /// Return the scalar product `x`ᵀ `y` of the vectors `x` and
    /// `y`.
    #[doc(alias = "gsl_blas_sdot")]
    pub fn dot(x: &VectorF32, y: &VectorF32) -> Result<f32, Error> {
        let mut result = 0.;
        let ret = unsafe { sys::gsl_blas_sdot(x.unwrap_shared(), y.unwrap_shared(), &mut result) };
        Error::handle(ret, result)
    }

    /// Return the scalar product `x`ᵀ `y` of the vectors `x` and
    /// `y`.
    #[doc(alias = "gsl_blas_dsdot")]
    pub fn ddot(x: &VectorF32, y: &VectorF32) -> Result<f64, Error> {
        let mut result = 0.;
        let ret = unsafe { sys::gsl_blas_dsdot(x.unwrap_shared(), y.unwrap_shared(), &mut result) };
        Error::handle(ret, result)
    }

    /// Return the Euclidean norm $‖x‖₂ = √{∑ x_i^2}$ of
    /// the vector `x`.
    #[doc(alias = "gsl_blas_snrm2")]
    pub fn nrm2(x: &VectorF32) -> f32 {
        unsafe { sys::gsl_blas_snrm2(x.unwrap_shared()) }
    }

    /// Return the absolute sum $∑ |x_i|$ of the elements of the
    /// vector `x`.
    #[doc(alias = "gsl_blas_sasum")]
    pub fn asum(x: &VectorF32) -> f32 {
        unsafe { sys::gsl_blas_sasum(x.unwrap_shared()) }
    }

    /// Return the index of the largest element of the vector `x`.
    /// The largest element is determined by its absolute
    /// magnitude .  If the largest value occurs several times
    /// then the index of the first occurrence is returned.
    #[doc(alias = "gsl_blas_isamax")]
    pub fn iamax(x: &VectorF32) -> usize {
        unsafe { sys::gsl_blas_isamax(x.unwrap_shared()) }
    }

    /// This function exchanges the elements of the vectors `x` and `y`.
    #[doc(alias = "gsl_blas_sswap")]
    pub fn swap(x: &mut VectorF32, y: &mut VectorF32) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_blas_sswap(x.unwrap_unique(), y.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// This function copy the elements of the vector `x` into the
    /// vector `y`.
    #[doc(alias = "gsl_blas_scopy")]
    pub fn copy(x: &mut VectorF32, y: &mut VectorF32) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_blas_scopy(x.unwrap_unique(), y.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// This function computes the sum `y` = `alpha` `x` + `y` for
    /// the vectors `x` and `y`.
    #[doc(alias = "gsl_blas_saxpy")]
    pub fn axpy(alpha: f32, x: &VectorF32, y: &mut VectorF32) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_blas_saxpy(alpha, x.unwrap_shared(), y.unwrap_unique()) };
        Error::handle(ret, ())
    }
    /// This function rescales the vector `x` by the
    /// multiplicative factor `alpha`.
    #[doc(alias = "gsl_blas_sscal")]
    pub fn scal(alpha: f32, x: &mut VectorF32) {
        unsafe { sys::gsl_blas_sscal(alpha, x.unwrap_unique()) }
    }

    /// This function computes a Givens rotation (c,s) which
    /// zeroes the vector (a,b),
    ///
    /// ```text
    /// [  c  s ] [ a ] = [ r ]
    /// [ -s  c ] [ b ]   [ 0 ]
    /// ```
    ///
    /// Return `(c, s, r)`.
    #[doc(alias = "gsl_blas_srotg")]
    pub fn rotg(mut a: f32, mut b: f32) -> Result<(f32, f32, f32), Error> {
        let mut c = 0.;
        let mut s = 0.;
        let ret = unsafe { sys::gsl_blas_srotg(&mut a, &mut b, &mut c, &mut s) };
        Error::handle(ret, (c, s, a))
    }
    /// This function applies a Givens rotation (x', y') = (c x + s y,
    /// -s x + c y) to the vectors x, y.
    #[doc(alias = "gsl_blas_srot")]
    pub fn rot(a: &mut VectorF32, b: &mut VectorF32, c: f32, d: f32) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_blas_srot(a.unwrap_unique(), b.unwrap_unique(), c, d) };
        Error::handle(ret, ())
    }

    /// Return a modified Givens transformation.  The modified Givens
    /// transformation is defined in the original [Level-1 BLAS
    /// specification](https://help.imsl.com/fortran/fnlmath/current/basic-linear-algebra-sub.htm#mch9_1817247609_srotmg).
    #[doc(alias = "gsl_blas_srotmg")]
    pub fn rotmg(mut d1: f32, mut d2: f32, mut b1: f32, b2: f32) -> Result<[f32; 5], Error> {
        let mut p = [f32::NAN; 5];
        let ret = unsafe { sys::gsl_blas_srotmg(&mut d1, &mut d2, &mut b1, b2, p.as_mut_ptr()) };
        Error::handle(ret, p)
    }

    /// This function applies a modified Givens transformation.
    #[doc(alias = "gsl_blas_srotm")]
    pub fn rotm(x: &mut VectorF32, y: &mut VectorF32, p: [f32; 5]) -> Result<(), Error> {
        let lenx = VectorF32::len(x);
        let leny = VectorF32::len(y);
        if lenx != leny {
            panic!("rgsl::blas::srotm: len(x) = {lenx} != len(y) = {leny}")
        }
        let ret = unsafe { sys::gsl_blas_srotm(x.unwrap_unique(), y.unwrap_unique(), p.as_ptr()) };
        Error::handle(ret, ())
    }

    // Level 2

    /// This function  computes the matrix-vector product and sum y
    /// = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for
    /// TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    #[doc(alias = "gsl_blas_sgemv")]
    pub fn gemv(
        transA: Transpose,
        alpha: f32,
        A: &MatrixF32,
        x: &VectorF32,
        beta: f32,
        y: &mut VectorF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_sgemv(
                transA.into(),
                alpha,
                A.unwrap_shared(),
                x.unwrap_shared(),
                beta,
                y.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-vector product x = op(A)
    /// x for the triang ular matrix A, where op(A) = A, A^T, A^H for
    /// TransA = CblasNoTrans, CblasTrans, CblasConjTrans.  When Uplo
    /// is CblasUpper then the upper triangle of A is used, and when
    /// Uplo is CblasLower then the lower triangle of A is used.  If
    /// Diag is CblasNonUnit then the diagonal of the matrix is used,
    /// but if Diag is CblasUnit then the diagonal elements of the
    /// matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_strmv")]
    pub fn trmv(
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        A: &MatrixF32,
        x: &mut VectorF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_strmv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes inv(op(A)) x for x, where op(A) = A,
    /// A^T, A^H for TransA = CblasNoTrans, CblasTrans,
    /// CblasConjTrans.  When Uplo is CblasUpper then the upper
    /// triangle of A is used, and when Uplo is CblasLower then the
    /// lower triangle of A is used.  If Diag is CblasNonUnit then the
    /// diagonal of the matrix is used, but if Diag is CblasUnit then
    /// the diagonal elements of the matrix A are taken as unity and
    /// are not referenced.
    #[doc(alias = "gsl_blas_strsv")]
    pub fn trsv(
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        A: &MatrixF32,
        x: &mut VectorF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_strsv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// These functions compute the matrix-vector product and sum y =
    /// \alpha A x + \beta y for the symmetric matrix A.  Since the
    /// matrix A is symmetric only its upper half or lower half need
    /// to be stored.  When Uplo is CblasUpper then the upper triangle
    /// and diagonal of A are used, and when Uplo is CblasLower then
    /// the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_ssymv")]
    pub fn symv(
        uplo: Uplo,
        alpha: f32,
        A: &MatrixF32,
        x: &VectorF32,
        beta: f32,
        y: &mut VectorF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ssymv(
                uplo.into(),
                alpha,
                A.unwrap_shared(),
                x.unwrap_shared(),
                beta,
                y.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A
    /// of the matrix A.
    #[doc(alias = "gsl_blas_sger")]
    pub fn ger(alpha: f32, x: &VectorF32, y: &VectorF32, A: &mut MatrixF32) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_sger(
                alpha,
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the symmetric rank-1 update A = \alpha
    /// x x^T + A of the symmetric matrix A. Since the matrix A is
    /// symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal
    /// of A are used, and when Uplo is CblasLower then the lower
    /// triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_ssyr")]
    pub fn syr(uplo: Uplo, alpha: f32, x: &VectorF32, A: &mut MatrixF32) -> Result<(), Error> {
        let ret =
            unsafe { sys::gsl_blas_ssyr(uplo.into(), alpha, x.unwrap_shared(), A.unwrap_unique()) };
        Error::handle(ret, ())
    }
    /// These functions compute the symmetric rank-2 update A = \alpha
    /// x y^T + \alpha y x^T + A of the symmetric matrix A.  Since the
    /// matrix A is symmetric only its upper half or lower half need
    /// to be stored.  When Uplo is CblasUpper then the upper triangle
    /// and diagonal of A are used, and when Uplo is CblasLower then
    /// the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_ssyr2")]
    pub fn syr2(
        uplo: Uplo,
        alpha: f32,
        x: &VectorF32,
        y: &VectorF32,
        A: &mut MatrixF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ssyr2(
                uplo.into(),
                alpha,
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }

    // Level 3

    /// This function computes the matrix-matrix product and sum C
    /// = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H
    /// for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and
    /// similarly for the parameter TransB.
    #[doc(alias = "gsl_blas_sgemm")]
    pub fn gemm(
        transA: Transpose,
        transB: Transpose,
        alpha: f32,
        A: &MatrixF32,
        B: &MatrixF32,
        beta: f32,
        C: &mut MatrixF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_sgemm(
                transA.into(),
                transB.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-matrix product and sum $C =
    /// α A B + β C$ for Side is CblasLeft and $C = α B A + β C$ for
    /// Side is CblasRight, where the matrix A is symmetric.  When
    /// Uplo is CblasUpper then the upper triangle and diagonal of A
    /// are used, and when Uplo is CblasLower then the lower triangle
    /// and diagonal of A are used.
    #[doc(alias = "gsl_blas_ssymm")]
    pub fn symm(
        side: Side,
        uplo: Uplo,
        alpha: f32,
        A: &MatrixF32,
        B: &MatrixF32,
        beta: f32,
        C: &mut MatrixF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ssymm(
                side.into(),
                uplo.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-matrix product B = \alpha
    /// op(A) B for Side is Left and B = \alpha B op(A) for Side is
    /// CblasRight.  The matrix A is triangular and op(A) = A, A^T,
    /// A^H for TransA = NoTrans, Trans, ConjTrans.  When Uplo is
    /// Upper then the upper triangle of A is used, and when Uplo is
    /// Lower then the lower triangle of A is used.  If Diag is
    /// NonUnit then the diagonal of A is used, but if Diag is Unit
    /// then the diagonal elements of the matrix A are taken as unity
    /// and are not referenced.
    #[doc(alias = "gsl_blas_strmm")]
    pub fn trmm(
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        alpha: f32,
        A: &MatrixF32,
        B: &mut MatrixF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_strmm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }

    /// This function computes the inverse-matrix matrix product B =
    /// \alpha op(inv(A))B for Side is Left and B = \alpha B
    /// op(inv(A)) for Side is Right.  The matrix A is triangular and
    /// op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and
    /// when Uplo is Lower then the lower triangle of A is used.  If
    /// Diag is NonUnit then the diagonal of A is used, but if Diag is
    /// Unit then the diagonal elements of the matrix A are taken as
    /// unity and are not referenced.
    #[doc(alias = "gsl_blas_strsm")]
    pub fn trsm(
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        alpha: f32,
        A: &MatrixF32,
        B: &mut MatrixF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_strsm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes a rank-k update of the symmetric matrix
    /// C, C = \alpha A A^T + \beta C when Trans is NoTrans and C =
    /// \alpha A^T A + \beta C when Trans is Trans.  Since the matrix
    /// C is symmetric only its upper half or lower half need to be
    /// stored.  When Uplo is Upper then the upper triangle and
    /// diagonal of C are used, and when Uplo is Lower then the lower
    /// triangle and diagonal of C are used.
    #[doc(alias = "gsl_blas_ssyrk")]
    pub fn syrk(
        uplo: Uplo,
        trans: Transpose,
        alpha: f32,
        A: &MatrixF32,
        beta: f32,
        C: &mut MatrixF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ssyrk(
                uplo.into(),
                trans.into(),
                alpha,
                A.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes a rank-2k update of the symmetric
    /// matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans
    /// is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when
    /// Trans is Trans.  Since the matrix C is symmetric only its
    /// upper half or lower half need to be stored.  When Uplo is
    /// Upper then the upper triangle and diagonal of C are used, and
    /// when Uplo is Lower then the lower triangle and diagonal of C
    /// are used.
    #[doc(alias = "gsl_blas_ssyr2k")]
    pub fn syr2k(
        uplo: Uplo,
        trans: Transpose,
        alpha: f32,
        A: &MatrixF32,
        B: &MatrixF32,
        beta: f32,
        C: &mut MatrixF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ssyr2k(
                uplo.into(),
                trans.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
}

/// `f64` vectors.
pub mod d {
    use super::*;
    use crate::{ffi::FFI, Error, MatrixF64, VectorF64};

    // Level 1

    /// Return the scalar product `x`ᵀ `y` of the vectors `x` and `y`.
    #[doc(alias = "gsl_blas_ddot")]
    pub fn dot(x: &VectorF64, y: &VectorF64) -> Result<f64, Error> {
        let mut result = 0.;
        let ret = unsafe { sys::gsl_blas_ddot(x.unwrap_shared(), y.unwrap_shared(), &mut result) };
        Error::handle(ret, result)
    }

    /// Return the Euclidean norm $‖x‖₂ = √{∑ x_i^2}$ of
    /// the vector `x`.
    #[doc(alias = "gsl_blas_dnrm2")]
    pub fn nrm2(x: &VectorF64) -> f64 {
        unsafe { sys::gsl_blas_dnrm2(x.unwrap_shared()) }
    }

    /// Return the absolute sum $∑ |x_i|$ of the elements of the
    /// vector `x`.
    #[doc(alias = "gsl_blas_dasum")]
    pub fn asum(x: &VectorF64) -> f64 {
        unsafe { sys::gsl_blas_dasum(x.unwrap_shared()) }
    }

    /// Return the index of the largest element of the vector `x`.
    /// The largest element is determined by its absolute
    /// magnitude.  If the largest value occurs several times then
    /// the index of the first occurrence is returned.
    #[doc(alias = "gsl_blas_idamax")]
    pub fn iamax(x: &VectorF64) -> usize {
        unsafe { sys::gsl_blas_idamax(x.unwrap_shared()) }
    }

    /// This function exchanges the elements of the vectors `x` and `y` .
    #[doc(alias = "gsl_blas_dswap")]
    pub fn swap(x: &mut VectorF64, y: &mut VectorF64) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_blas_dswap(x.unwrap_unique(), y.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// This function copy the elements of the vector `x` into the
    /// vector `y`.
    #[doc(alias = "gsl_blas_dcopy")]
    pub fn copy(x: &mut VectorF64, y: &mut VectorF64) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_blas_dcopy(x.unwrap_unique(), y.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// This function computes the sum `y` = `alpha` `x` + `y` for
    /// the vectors `x` and `y`.
    #[doc(alias = "gsl_blas_daxpy")]
    pub fn axpy(alpha: f64, x: &VectorF64, y: &mut VectorF64) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_blas_daxpy(alpha, x.unwrap_shared(), y.unwrap_unique()) };
        Error::handle(ret, ())
    }
    /// This function rescales the vector `x` by the
    /// multiplicative factor `alpha`.
    #[doc(alias = "gsl_blas_dscal")]
    pub fn scal(alpha: f64, x: &mut VectorF64) {
        unsafe { sys::gsl_blas_dscal(alpha, x.unwrap_unique()) }
    }

    /// This function computes a Givens rotation (c,s) which
    /// zeroes the vector (a,b),
    ///
    /// ```text
    /// [  c  s ] [ a ] = [ r ]
    /// [ -s  c ] [ b ]   [ 0 ]
    /// ```
    ///
    /// Return `(c, s, r)`.
    #[doc(alias = "gsl_blas_drotg")]
    pub fn rotg(mut a: f64, mut b: f64) -> Result<(f64, f64, f64), Error> {
        let mut c = 0.;
        let mut s = 0.;
        let ret = unsafe { sys::gsl_blas_drotg(&mut a, &mut b, &mut c, &mut s) };
        Error::handle(ret, (c, s, a))
    }

    /// This function applies a Givens rotation (x', y') = (c x +
    /// s y, -s x + c y) to the vectors x, y.
    #[doc(alias = "gsl_blas_drot")]
    pub fn rot(a: &mut VectorF64, b: &mut VectorF64, c: f64, d: f64) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_blas_drot(a.unwrap_unique(), b.unwrap_unique(), c, d) };
        Error::handle(ret, ())
    }

    /// Return a modified Givens transformation.  The modified
    /// Givens transformation is defined in the original [Level-1
    /// BLAS
    /// specification](https://help.imsl.com/fortran/fnlmath/current/basic-linear-algebra-sub.htm#mch9_1817247609_srotmg).
    #[doc(alias = "gsl_blas_drotmg")]
    pub fn rotmg(mut d1: f64, mut d2: f64, mut b1: f64, b2: f64) -> Result<[f64; 5], Error> {
        let mut p = [f64::NAN; 5];
        let ret = unsafe { sys::gsl_blas_drotmg(&mut d1, &mut d2, &mut b1, b2, p.as_mut_ptr()) };
        Error::handle(ret, p)
    }

    /// This function applies a modified Givens transformation.
    #[doc(alias = "gsl_blas_drotm")]
    pub fn drotm(x: &mut VectorF64, y: &mut VectorF64, p: [f64; 5]) -> Result<(), Error> {
        let lenx = VectorF64::len(x);
        let leny = VectorF64::len(y);
        if lenx != leny {
            panic!("rgsl::blas::drotm: len(x) = {lenx} != len(y) = {leny}")
        }
        let ret = unsafe { sys::gsl_blas_drotm(x.unwrap_unique(), y.unwrap_unique(), p.as_ptr()) };
        Error::handle(ret, ())
    }

    // Level 2

    /// This function computes the matrix-vector product and sum y
    /// = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for
    /// TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    #[doc(alias = "gsl_blas_dgemv")]
    pub fn gemv(
        transA: Transpose,
        alpha: f64,
        A: &MatrixF64,
        x: &VectorF64,
        beta: f64,
        y: &mut VectorF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_dgemv(
                transA.into(),
                alpha,
                A.unwrap_shared(),
                x.unwrap_shared(),
                beta,
                y.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-vector product x = op(A) x
    /// for the triangular matrix A, where op(A) = A, A^T, A^H for
    /// TransA = CblasNoTrans, CblasTrans, CblasConjTrans.  When Uplo
    /// is CblasUpper then the upper triangle of A is used, and when
    /// Uplo is CblasLower then the lower triangle of A is used.  If
    /// Diag is CblasNonUnit then the diagonal of the matrix is used,
    /// but if Diag is CblasUnit then the diagonal elements of the
    /// matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_dtrmv")]
    pub fn trmv(
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        A: &MatrixF64,
        x: &mut VectorF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_dtrmv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes inv(op(A)) x for x, where op(A) = A,
    /// A^T, A^H for TransA = CblasNoTrans, CblasTrans,
    /// CblasConjTrans.  When Uplo is CblasUpper then the upper
    /// triangle of A is used, and when Uplo is CblasLower then the
    /// lower triangle of A is used.  If Diag is CblasNonUnit then the
    /// diagonal of the matrix is used, but if Diag is CblasUnit then
    /// the diagonal elements of the matrix A are taken as unity and
    /// are not referenced.
    #[doc(alias = "gsl_blas_dtrsv")]
    pub fn trsv(
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        A: &MatrixF64,
        x: &mut VectorF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_dtrsv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// These functions compute the matrix-vector product and sum y =
    /// \alpha A x + \beta y for the symmetric matrix A.  Since the
    /// matrix A is symmetric only its upper half or lower half need
    /// to be stored.  When Uplo is CblasUpper then the upper triangle
    /// and diagonal of A are used, and when Uplo is CblasLower then
    /// the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_dsymv")]
    pub fn symv(
        uplo: Uplo,
        alpha: f64,
        A: &MatrixF64,
        x: &VectorF64,
        beta: f64,
        y: &mut VectorF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_dsymv(
                uplo.into(),
                alpha,
                A.unwrap_shared(),
                x.unwrap_shared(),
                beta,
                y.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the rank-1 update A = \alpha x y^T + A
    /// of the matrix A.
    #[doc(alias = "gsl_blas_dger")]
    pub fn ger(alpha: f64, x: &VectorF64, y: &VectorF64, A: &mut MatrixF64) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_dger(
                alpha,
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the symmetric rank-1 update A = \alpha
    /// x x^T + A of the symmetric matrix A. Since the matrix A is
    /// symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal
    /// of A are used, and when Uplo is CblasLower then the lower
    /// triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_dsyr")]
    pub fn syr(uplo: Uplo, alpha: f64, x: &VectorF64, A: &mut MatrixF64) -> Result<(), Error> {
        let ret =
            unsafe { sys::gsl_blas_dsyr(uplo.into(), alpha, x.unwrap_shared(), A.unwrap_unique()) };
        Error::handle(ret, ())
    }
    /// These functions compute the symmetric rank-2 update A = \alpha
    /// x y^T + \alpha y x^T + A of the symmetric matrix A.  Since the
    /// matrix A is symmetric only its upper half or lower half need
    /// to be stored.  When Uplo is CblasUpper then the upper triangle
    /// and diagonal of A are used, and when Uplo is CblasLower then
    /// the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_dsyr2")]
    pub fn syr2(
        uplo: Uplo,
        alpha: f64,
        x: &VectorF64,
        y: &VectorF64,
        A: &mut MatrixF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_dsyr2(
                uplo.into(),
                alpha,
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }

    // Level 3

    /// This function computes the matrix-matrix product and sum C
    /// = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H
    /// for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and
    /// similarly for the parameter TransB.
    #[doc(alias = "gsl_blas_dgemm")]
    pub fn gemm(
        transA: Transpose,
        transB: Transpose,
        alpha: f64,
        A: &MatrixF64,
        B: &MatrixF64,
        beta: f64,
        C: &mut MatrixF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_dgemm(
                transA.into(),
                transB.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-matrix product and sum $C =
    /// α A B + β C$ for Side is `CblasLeft` and $C = α B A + β C$ for
    /// Side is `CblasRight`, where the matrix $A$ is symmetric.  When
    /// Uplo is CblasUpper then the upper triangle and diagonal of A
    /// are used, and when Uplo is CblasLower then the lower triangle
    /// and diagonal of A are used.
    #[doc(alias = "gsl_blas_dsymm")]
    pub fn symm(
        side: Side,
        uplo: Uplo,
        alpha: f64,
        A: &MatrixF64,
        B: &MatrixF64,
        beta: f64,
        C: &mut MatrixF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_dsymm(
                side.into(),
                uplo.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-matrix product B = \alpha
    /// op(A) B for Side is Left and B = \alpha B op(A) for Side is
    /// CblasRight.  The matrix A is triangular and op(A) = A, A^T,
    /// A^H for TransA = NoTrans, Trans, ConjTrans.  When Uplo is
    /// Upper then the upper triangle of A is used, and when Uplo is
    /// Lower then the lower triangle of A is used.  If Diag is
    /// NonUnit then the diagonal of A is used, but if Diag is Unit
    /// then the diagonal elements of the matrix A are taken as unity
    /// and are not referenced.
    #[doc(alias = "gsl_blas_dtrmm")]
    pub fn ddtrmm(
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        alpha: f64,
        A: &MatrixF64,
        B: &mut MatrixF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_dtrmm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the inverse-matrix matrix product B =
    /// \alpha op(inv(A))B for Side is Left and B = \alpha B
    /// op(inv(A)) for Side is Right.  The matrix A is triangular and
    /// op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and
    /// when Uplo is Lower then the lower triangle of A is used.  If
    /// Diag is NonUnit then the diagonal of A is used, but if Diag is
    /// Unit then the di agonal elements of the matrix A are taken as
    /// unity and are not referenced.
    #[doc(alias = "gsl_blas_dtrsm")]
    pub fn trsm(
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        alpha: f64,
        A: &MatrixF64,
        B: &mut MatrixF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_dtrsm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes a rank-k update of the symmetric matrix
    /// C, C = \alpha A A^T + \beta C when Trans is NoTrans and C =
    /// \alpha A^T A + \beta C when Trans is Trans.  Since the matrix
    /// C is symmetric only its upper half or lower half need to be
    /// stored.  When Uplo is Upper then the upper triangle and
    /// diagonal of C are used, and when Uplo is Lower then the lower
    /// triangle and diagonal of C are used.
    #[doc(alias = "gsl_blas_dsyrk")]
    pub fn syrk(
        uplo: Uplo,
        trans: Transpose,
        alpha: f64,
        A: &MatrixF64,
        beta: f64,
        C: &mut MatrixF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_dsyrk(
                uplo.into(),
                trans.into(),
                alpha,
                A.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes a rank-2k update of the symmetric
    /// matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans
    /// is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when
    /// Trans is Trans.  Since the matrix C is symmetric only its
    /// upper half or lower half need to be stored.  When Uplo is
    /// Upper then the upper triangle and diagonal of C are used, and
    /// when Uplo is Lower then the lower triangle and diagonal of C
    /// are used.
    #[doc(alias = "gsl_blas_dsyr2k")]
    pub fn syr2k(
        uplo: Uplo,
        trans: Transpose,
        alpha: f64,
        A: &MatrixF64,
        B: &MatrixF64,
        beta: f64,
        C: &mut MatrixF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_dsyr2k(
                uplo.into(),
                trans.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
}

/// `Complex<f32>` vectors.
#[cfg(feature = "complex")]
#[cfg_attr(docsrs, doc(cfg(feature = "complex")))]
pub mod c {
    use super::*;
    use crate::{
        ffi::FFI,
        types::complex::{FromC, ToC},
    };
    use crate::{Error, MatrixComplexF32, VectorComplexF32};
    use num_complex::Complex;

    // Level 1

    /// Return the complex scalar product `x`ᵀ y for the
    /// vectors `x` and `y`.
    #[doc(alias = "gsl_blas_cdotu")]
    pub fn dotu(x: &VectorComplexF32, y: &VectorComplexF32) -> Result<Complex<f32>, Error> {
        let mut dotu = Complex::<f32>::default().unwrap();
        let ret = unsafe { sys::gsl_blas_cdotu(x.unwrap_shared(), y.unwrap_shared(), &mut dotu) };
        Error::handle(ret, dotu.wrap())
    }

    /// Return the complex conjugate scalar product `x`ᴴ `y` for
    /// the vectors `x` and `y`.
    #[doc(alias = "gsl_blas_cdotc")]
    pub fn dot(x: &VectorComplexF32, y: &VectorComplexF32) -> Result<Complex<f32>, Error> {
        let mut dotc = Complex::<f32>::default().unwrap();
        let ret = unsafe { sys::gsl_blas_cdotc(x.unwrap_shared(), y.unwrap_shared(), &mut dotc) };
        Error::handle(ret, dotc.wrap())
    }

    /// Return the Euclidean norm of the complex vector `x`,
    /// $‖x‖_2 = √{∑ (\Re(x_i)^2 + \Im(x_i)^2)}$.
    #[doc(alias = "gsl_blas_scnrm2")]
    pub fn nrm2(x: &VectorComplexF32) -> f32 {
        unsafe { sys::gsl_blas_scnrm2(x.unwrap_shared()) }
    }

    /// Return the sum of the magnitudes of the real and imaginary
    /// parts of the complex vector `x`, $∑ |\Re(x_i)| +
    /// |\Im(x_i)|$.
    #[doc(alias = "gsl_blas_scasum")]
    pub fn asum(x: &VectorComplexF32) -> f32 {
        unsafe { sys::gsl_blas_scasum(x.unwrap_shared()) }
    }

    /// Return the index of the largest element of the vector `x`.
    /// The largest element is determined by the sum of the
    /// magnitudes of the real and imaginary parts $|\Re(x_i)| +
    /// |\Im(x_i)|$ for complex vectors.  If the largest value
    /// occurs several times then the index of the first
    /// occurrence is returned.
    #[doc(alias = "gsl_blas_icamax")]
    pub fn iamax(x: &VectorComplexF32) -> usize {
        unsafe { sys::gsl_blas_icamax(x.unwrap_shared()) }
    }

    /// This function exchanges the elements of the vectors `x` and `y`.
    #[doc(alias = "gsl_blas_cswap")]
    pub fn cswap(x: &mut VectorComplexF32, y: &mut VectorComplexF32) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_blas_cswap(x.unwrap_unique(), y.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// This function copy the elements of the vector `x` into the
    /// vector `y`.
    #[doc(alias = "gsl_blas_ccopy")]
    pub fn copy(x: &mut VectorComplexF32, y: &mut VectorComplexF32) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_blas_ccopy(x.unwrap_unique(), y.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// This function computes the sum `y` = `alpha` `x` + `y` for
    /// the vectors `x` and `y`.
    #[doc(alias = "gsl_blas_caxpy")]
    pub fn caxpy(
        alpha: &Complex<f32>,
        x: &VectorComplexF32,
        y: &mut VectorComplexF32,
    ) -> Result<(), Error> {
        let ret =
            unsafe { sys::gsl_blas_caxpy(alpha.unwrap(), x.unwrap_shared(), y.unwrap_unique()) };
        Error::handle(ret, ())
    }
    /// This function rescales the vector `x` by the multiplicative
    /// factor `alpha`.
    #[doc(alias = "gsl_blas_cscal")]
    pub fn scal(alpha: &Complex<f32>, x: &mut VectorComplexF32) {
        unsafe { sys::gsl_blas_cscal(alpha.unwrap(), x.unwrap_unique()) }
    }

    /// This function rescales the vector `x` by the multiplicative
    /// factor `alpha`.
    #[doc(alias = "gsl_blas_csscal")]
    pub fn rscal(alpha: f32, x: &mut VectorComplexF32) {
        unsafe { sys::gsl_blas_csscal(alpha, x.unwrap_unique()) }
    }

    // Level 2

    /// This function computes the matrix-vector product and sum y
    /// = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for
    /// TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    #[doc(alias = "gsl_blas_cgemv")]
    pub fn gemv(
        transA: Transpose,
        alpha: Complex<f32>,
        A: &MatrixComplexF32,
        x: &VectorComplexF32,
        beta: &Complex<f32>,
        y: &mut VectorComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_cgemv(
                transA.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                x.unwrap_shared(),
                beta.unwrap(),
                y.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-vector product x = op(A) x
    /// for the triangular matrix A, where op(A) = A, A^T, A^H for
    /// TransA = CblasNoTrans, CblasTrans, CblasConjTrans.  When Uplo
    /// is CblasUpper then the upper triangle of A is used, and when
    /// Uplo is CblasLower then the lower triangle of A is used.  If
    /// Diag is CblasNonUnit then the diagonal of the matrix is used,
    /// but if Diag is CblasUnit then the diagonal elements of the
    /// matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_ctrmv")]
    pub fn trmv(
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        A: &MatrixComplexF32,
        x: &mut VectorComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ctrmv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes inv(op(A)) x for x, where op(A) = A,
    /// A^T, A^H for TransA = CblasNoTrans, CblasTrans,
    /// CblasConjTrans.  When Uplo is CblasUpper then the upper
    /// triangle of A is used, and when Uplo is CblasLower then the
    /// lower triangle of A is used.  If Diag is CblasNonUnit then the
    /// diagonal of the matrix is used, but if Diag is CblasUnit then
    /// the diagonal elements of the matrix A are taken as unity and
    /// are not referenced.
    #[doc(alias = "gsl_blas_ctrsv")]
    pub fn trsv(
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        A: &MatrixComplexF32,
        x: &mut VectorComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ctrsv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// These functions compute the matrix-vector product and sum y =
    /// \alpha A x + \beta y for the hermitian matrix A.  Since the
    /// matrix A is hermitian only its upper half or lower half need
    /// to be stored. When Uplo is CblasUpper then the upper triangle
    /// and diagonal of A are used, and when Uplo is CblasLower then
    /// the lower triangle and diagonal of A are used.  The imaginary
    /// elements of the diagonal are automatically assumed to be zero
    /// and are not referenced.
    #[doc(alias = "gsl_blas_chemv")]
    pub fn hemv(
        uplo: Uplo,
        alpha: &Complex<f32>,
        A: &MatrixComplexF32,
        x: &VectorComplexF32,
        beta: &Complex<f32>,
        y: &mut VectorComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_chemv(
                uplo.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                x.unwrap_shared(),
                beta.unwrap(),
                y.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the rank-1 update A = \alpha x y^T + A
    /// of the matrix A.
    #[doc(alias = "gsl_blas_cgeru")]
    pub fn geru(
        alpha: &Complex<f32>,
        x: &VectorComplexF32,
        y: &VectorComplexF32,
        A: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_cgeru(
                alpha.unwrap(),
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the conjugate rank-1 update A = \alpha
    /// x y^H + A of the matrix A.
    #[doc(alias = "gsl_blas_cgerc")]
    pub fn gerc(
        alpha: &Complex<f32>,
        x: &VectorComplexF32,
        y: &VectorComplexF32,
        A: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_cgerc(
                alpha.unwrap(),
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }

    /// These functions compute the hermitian rank-1 update A = \alpha
    /// x x^H + A of the hermitian matrix A.  Since the matrix A is
    /// hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal
    /// of A are used, and when Uplo is CblasLower then the lower
    /// triangle and diagonal of A are used.  The imaginary elements
    /// of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_cher")]
    pub fn her(
        uplo: Uplo,
        alpha: f32,
        x: &VectorComplexF32,
        A: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret =
            unsafe { sys::gsl_blas_cher(uplo.into(), alpha, x.unwrap_shared(), A.unwrap_unique()) };
        Error::handle(ret, ())
    }
    /// These functions compute the hermitian rank-2 update A = \alpha
    /// x y^H + \alpha^* y x^H + A of the hermitian matrix A.  Since
    /// the matrix A is hermitian only its upper half or lower half
    /// need to be stored.  When Uplo is CblasUpper then the upper
    /// triangle and diagonal of A are used, and when Uplo is
    /// CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set
    /// to zero.
    #[doc(alias = "gsl_blas_cher2")]
    pub fn her2(
        uplo: Uplo,
        alpha: &Complex<f32>,
        x: &VectorComplexF32,
        y: &VectorComplexF32,
        A: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_cher2(
                uplo.into(),
                alpha.unwrap(),
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-matrix product and sum C
    /// = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H
    /// for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and
    /// similarly for the parameter TransB.
    #[doc(alias = "gsl_blas_cgemm")]
    pub fn gemm(
        transA: Transpose,
        transB: Transpose,
        alpha: &Complex<f32>,
        A: &MatrixComplexF32,
        B: &MatrixComplexF32,
        beta: &Complex<f32>,
        C: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_cgemm(
                transA.into(),
                transB.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta.unwrap(),
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-matrix product and sum $C =
    /// α A B + β C$ for Side is `CblasLeft` and $C = α B A + beta C$
    /// for Side is `CblasRight`, where the matrix $A$ is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal
    /// of A are used, and when Uplo is CblasLower then the lower
    /// triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_csymm")]
    pub fn symm(
        side: Side,
        uplo: Uplo,
        alpha: &Complex<f32>,
        A: &MatrixComplexF32,
        B: &MatrixComplexF32,
        beta: &Complex<f32>,
        C: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_csymm(
                side.into(),
                uplo.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta.unwrap(),
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-matrix product and sum C =
    /// \alpha A B + \beta C for Side is Left and C = \alpha B A +
    /// \beta C for Side is Right, where the matrix A is hermitian.
    /// When Uplo is Upper then the upper triangle and diagonal of A
    /// are used, and when Uplo is Lower then the lower triangle and
    /// diagonal of A are used.  The imaginary elements of the
    /// diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_chemm")]
    pub fn hemm(
        side: Side,
        uplo: Uplo,
        alpha: &Complex<f32>,
        A: &MatrixComplexF32,
        B: &MatrixComplexF32,
        beta: &Complex<f32>,
        C: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_chemm(
                side.into(),
                uplo.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta.unwrap(),
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-matrix product B = \alpha
    /// op(A) B for Side is Left and B = \alpha B op(A) for Side is
    /// CblasRight.  The matrix A is triangular and op(A) = A, A^T,
    /// A^H for TransA = NoTrans, Trans, ConjTrans.  When Uplo is
    /// Upper then the upper triangle of A is used, and when Uplo is
    /// Lower then the lower triangle of A is used.  If Diag is
    /// NonUnit then the diagonal of A is used, but if Diag is Unit
    /// then the diagonal elements of the matrix A are taken as unity
    /// and are not referenced.
    #[doc(alias = "gsl_blas_ctrmm")]
    pub fn trmm(
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        alpha: &Complex<f32>,
        A: &MatrixComplexF32,
        B: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ctrmm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the inverse-matrix matrix product B =
    /// \alpha op(inv(A))B for Side is Left and B = \alpha B
    /// op(inv(A)) for Side is Right.  The matrix A is triangular and
    /// op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and
    /// when Uplo is Lower then the lower triangle of A is used.  If
    /// Diag is NonUnit then the diagonal of A is used, but if Diag is
    /// Unit then the diagonal elements of the matrix A are taken as
    /// unity and are not referenced.
    #[doc(alias = "gsl_blas_ctrsm")]
    pub fn trsm(
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        alpha: &Complex<f32>,
        A: &MatrixComplexF32,
        B: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ctrsm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes a rank-k update of the symmetric matrix
    /// C, C = \alpha A A^T + \beta C when Trans is NoTrans and C =
    /// \alpha A^T A + \beta C when Trans is Trans.  Since the matrix
    /// C is symmetric only its upper half or lower half need to be
    /// stored.  When Uplo is Upper then the upper triangle and
    /// diagonal of C are used, and when Uplo is Lower then the lower
    /// triangle and diagonal of C are used.
    #[doc(alias = "gsl_blas_csyrk")]
    pub fn syrk(
        uplo: Uplo,
        trans: Transpose,
        alpha: &Complex<f32>,
        A: &MatrixComplexF32,
        beta: &Complex<f32>,
        C: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_csyrk(
                uplo.into(),
                trans.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                beta.unwrap(),
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// These functions compute a rank-k update of the hermitian
    /// matrix C, C = \alpha A A^H + \beta C when Trans is NoTrans and
    /// C = \alpha A^H A + \beta C when Trans is ConjTrans.  Since the
    /// matrix C is hermitian only its upper half or lower half need
    /// to be stored.  When Uplo is Upper then the upper triangle and
    /// diagonal of C are used, and when Uplo is Lower then the lower
    /// triangle and diagonal of C are used.  The imaginary elements
    /// of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_cherk")]
    pub fn herk(
        uplo: Uplo,
        trans: Transpose,
        alpha: f32,
        A: &MatrixComplexF32,
        beta: f32,
        C: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_cherk(
                uplo.into(),
                trans.into(),
                alpha,
                A.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes a rank-2k update of the symmetric
    /// matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans
    /// is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when
    /// Trans is Trans.  Since the matrix C is symmetric only its
    /// upper half or lower half need to be stored.  When Uplo is
    /// Upper then the upper triangle and diagonal of C are used, and
    /// when Uplo is Lower then the lower triangle and diagonal of C
    /// are used.
    #[doc(alias = "gsl_blas_csyr2k")]
    #[cfg(feature = "complex")]
    pub fn syr2k(
        uplo: Uplo,
        trans: Transpose,
        alpha: &Complex<f32>,
        A: &MatrixComplexF32,
        B: &MatrixComplexF32,
        beta: &Complex<f32>,
        C: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_csyr2k(
                uplo.into(),
                trans.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta.unwrap(),
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes a rank-2k update of the hermitian
    /// matrix C, C = \alpha A B^H + \alpha^* B A^H + \beta C when
    /// Trans is NoTrans and C = \alpha A^H B + \alpha^* B^H A + \beta
    /// C when Trans is ConjTrans.  Since the matrix C is hermitian
    /// only its upper half or lower half need to be stored.  When
    /// Uplo is Upper then the upper triangle and diagonal of C are
    /// used, and when Uplo is Lower then the lower triangle and
    /// diagonal of C are used.  The imaginary elements of the
    /// diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_cher2k")]
    pub fn her2k(
        uplo: Uplo,
        trans: Transpose,
        alpha: &Complex<f32>,
        A: &MatrixComplexF32,
        B: &MatrixComplexF32,
        beta: f32,
        C: &mut MatrixComplexF32,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_cher2k(
                uplo.into(),
                trans.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
}

/// `Complex<f64>` vectors.
#[cfg(feature = "complex")]
#[cfg_attr(docsrs, doc(cfg(feature = "complex")))]
pub mod z {
    use super::*;
    use crate::{
        ffi::FFI,
        types::complex::{FromC, ToC},
    };
    use crate::{Error, MatrixComplexF64, VectorComplexF64};
    use num_complex::Complex;

    // Level 1

    /// Return the complex scalar product `x`ᵀ `y` for the
    /// vectors `x` and `y`.
    #[doc(alias = "gsl_blas_zdotu")]
    pub fn dotu(x: &VectorComplexF64, y: &VectorComplexF64) -> Result<Complex<f64>, Error> {
        let mut dotu = Complex::<f64>::default().unwrap();
        let ret = unsafe { sys::gsl_blas_zdotu(x.unwrap_shared(), y.unwrap_shared(), &mut dotu) };
        Error::handle(ret, dotu.wrap())
    }

    /// Return the complex conjugate scalar product `x`ᴴ `y` for
    /// the vectors `x` and `y`.
    #[doc(alias = "gsl_blas_zdotc")]
    pub fn dot(x: &VectorComplexF64, y: &VectorComplexF64) -> Result<Complex<f64>, Error> {
        let mut dotc = Complex::<f64>::default().unwrap();
        let ret = unsafe { sys::gsl_blas_zdotc(x.unwrap_shared(), y.unwrap_shared(), &mut dotc) };
        Error::handle(ret, dotc.wrap())
    }

    /// Return the Euclidean norm of the complex vector `x`,
    /// $‖x‖_2 = √{∑ (\Re(x_i)^2 + \Im(x_i)^2)}$.
    #[doc(alias = "gsl_blas_dznrm2")]
    pub fn nrm2(x: &VectorComplexF64) -> f64 {
        unsafe { sys::gsl_blas_dznrm2(x.unwrap_shared()) }
    }

    /// Return the sum of the magnitudes of the real and imaginary
    /// parts of the complex vector `x`, $∑ |\Re(x_i)| +
    /// |\Im(x_i)|$.
    #[doc(alias = "gsl_blas_dzasum")]
    pub fn asum(x: &VectorComplexF64) -> f64 {
        unsafe { sys::gsl_blas_dzasum(x.unwrap_shared()) }
    }

    /// Return the index of the largest element of the vector `x`.
    /// The largest element is determined by the sum of the
    /// magnitudes of the real and imaginary parts $|\Re(x_i)| +
    /// |\Im(x_i)|$ for complex vectors.  If the largest value
    /// occurs several times then the index of the first
    /// occurrence is returned.
    #[doc(alias = "gsl_blas_izamax")]
    pub fn iamax(x: &VectorComplexF64) -> usize {
        unsafe { sys::gsl_blas_izamax(x.unwrap_shared()) }
    }

    /// This function exchanges the elements of the vectors `x` and `y`.
    #[doc(alias = "gsl_blas_zswap")]
    pub fn swap(x: &mut VectorComplexF64, y: &mut VectorComplexF64) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_blas_zswap(x.unwrap_unique(), y.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// This function copy the elements of the vector `x` into the
    /// vector `y`.
    #[doc(alias = "gsl_blas_zcopy")]
    pub fn zcopy(x: &mut VectorComplexF64, y: &mut VectorComplexF64) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_blas_zcopy(x.unwrap_unique(), y.unwrap_unique()) };
        Error::handle(ret, ())
    }
    /// This function computes the sum `y` = `alpha` `x` + `y` for
    /// the vectors `x` and `y`.
    #[doc(alias = "gsl_blas_zaxpy")]
    pub fn axpy(
        alpha: &Complex<f64>,
        x: &VectorComplexF64,
        y: &mut VectorComplexF64,
    ) -> Result<(), Error> {
        let ret =
            unsafe { sys::gsl_blas_zaxpy(alpha.unwrap(), x.unwrap_shared(), y.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// This function rescales the vector `x` by the multiplicative
    /// factor `alpha`.
    #[doc(alias = "gsl_blas_zscal")]
    pub fn scal(alpha: &Complex<f64>, x: &mut VectorComplexF64) {
        unsafe { sys::gsl_blas_zscal(alpha.unwrap(), x.unwrap_unique()) }
    }

    /// This function rescales the vector `x` by the multiplicative
    /// factor `alpha`.
    #[doc(alias = "gsl_blas_zdscal")]
    pub fn rscal(alpha: f64, x: &mut VectorComplexF64) {
        unsafe { sys::gsl_blas_zdscal(alpha, x.unwrap_unique()) }
    }

    // Level 2

    /// This function computes the matrix-vector product and sum y
    /// = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for
    /// TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    #[doc(alias = "gsl_blas_zgemv")]
    pub fn gemv(
        transA: Transpose,
        alpha: &Complex<f64>,
        A: &MatrixComplexF64,
        x: &VectorComplexF64,
        beta: &Complex<f64>,
        y: &mut VectorComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_zgemv(
                transA.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                x.unwrap_shared(),
                beta.unwrap(),
                y.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-vector product x = op(A) x
    /// for the triangular matrix A, where op(A) = A, A^T, A^H for
    /// TransA = CblasNoTrans, CblasTrans, CblasConjTrans.  When Uplo
    /// is CblasUpper then the upper triangle of A is used, and when
    /// Uplo is CblasLower then the lower triangle of A is used.  If
    /// Diag is CblasNonUnit then the diagonal of the matrix is used,
    /// but if Diag is CblasUnit then the diagonal elements of the
    /// matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_ztrmv")]
    pub fn trmv(
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        A: &MatrixComplexF64,
        x: &mut VectorComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ztrmv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes inv(op(A)) x for x, where op(A) = A,
    /// A^T, A^H for TransA = CblasNoTrans, CblasTrans,
    /// CblasConjTrans.  When Uplo is CblasUpper then the upper
    /// triangle of A is used, and when Uplo is CblasLower then the
    /// lower triangle of A is used.  If Diag is CblasNonUnit then the
    /// diagonal of the matrix is used, but if Diag is CblasUnit then
    /// the diagonal elements of the matrix A are taken as unity and
    /// are not referenced.
    #[doc(alias = "gsl_blas_ztrsv")]
    pub fn trsv(
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        A: &MatrixComplexF64,
        x: &mut VectorComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ztrsv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// These functions compute the matrix-vector product and sum y =
    /// \alpha A x + \beta y for the hermitian matrix A.  Since the
    /// matrix A is hermitian only its upper half or lower half need
    /// to be stored. When Uplo is CblasUpper then the upper triangle
    /// and diagonal of A are used, and when Uplo is CblasLower then
    /// the lower triangle and diagonal of A are used.  The imaginary
    /// elements of the diagonal are automatically assumed to be zero
    /// and are not referenced.
    #[doc(alias = "gsl_blas_zhemv")]
    pub fn hemv(
        uplo: Uplo,
        alpha: &Complex<f64>,
        A: &MatrixComplexF64,
        x: &VectorComplexF64,
        beta: &Complex<f64>,
        y: &mut VectorComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_zhemv(
                uplo.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                x.unwrap_shared(),
                beta.unwrap(),
                y.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the rank-1 update A = \alpha x y^T + A
    /// of the matrix A.
    #[doc(alias = "gsl_blas_zgeru")]
    pub fn geru(
        alpha: &Complex<f64>,
        x: &VectorComplexF64,
        y: &VectorComplexF64,
        A: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_zgeru(
                alpha.unwrap(),
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the conjugate rank-1 update A = \alpha
    /// x y^H + A of the matrix A.
    #[doc(alias = "gsl_blas_zgerc")]
    pub fn gerc(
        alpha: &Complex<f64>,
        x: &VectorComplexF64,
        y: &VectorComplexF64,
        A: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_zgerc(
                alpha.unwrap(),
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// These functions compute the hermitian rank-1 update A = \alpha
    /// x x^H + A of the hermitian matrix A.  Since the matrix A is
    /// hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal
    /// of A are used, and when Uplo is CblasLower then the lower
    /// triangle and diagonal of A are used.  The imaginary elements
    /// of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_zher")]
    pub fn her(
        uplo: Uplo,
        alpha: f64,
        x: &VectorComplexF64,
        A: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret =
            unsafe { sys::gsl_blas_zher(uplo.into(), alpha, x.unwrap_shared(), A.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// These functions compute the hermitian rank-2 update A = \alpha
    /// x y^H + \alpha^* y x^H + A of the hermitian matrix A.  Since
    /// the matrix A is hermitian only its upper half or lower half
    /// need to be stored.  When Uplo is CblasUpper then the upper
    /// triangle and diagonal of A are used, and when Uplo is
    /// CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set
    /// to zero.
    #[doc(alias = "gsl_blas_zher2")]
    pub fn dzher2(
        uplo: Uplo,
        alpha: &Complex<f64>,
        x: &VectorComplexF64,
        y: &VectorComplexF64,
        A: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_zher2(
                uplo.into(),
                alpha.unwrap(),
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }

    // Level 3

    /// This function computes the matrix-matrix product and sum C
    /// = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H
    /// for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and
    /// similarly for the parameter TransB.
    #[doc(alias = "gsl_blas_zgemm")]
    pub fn gemm(
        transA: Transpose,
        transB: Transpose,
        alpha: &Complex<f64>,
        A: &MatrixComplexF64,
        B: &MatrixComplexF64,
        beta: &Complex<f64>,
        C: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_zgemm(
                transA.into(),
                transB.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta.unwrap(),
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-matrix product and sum $C =
    /// α A B + β C$ for Side is `CblasLeft` and `C = α B A + β C$ for
    /// Side is `CblasRight`, where the matrix $A$ is symmetric.  When
    /// Uplo is CblasUpper then the upper triangle and diagonal of A
    /// are used, and when Uplo is CblasLower then the lower triangle
    /// and diagonal of A are used.
    #[doc(alias = "gsl_blas_zsymm")]
    pub fn symm(
        side: Side,
        uplo: Uplo,
        alpha: &Complex<f64>,
        A: &MatrixComplexF64,
        B: &MatrixComplexF64,
        beta: &Complex<f64>,
        C: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_zsymm(
                side.into(),
                uplo.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta.unwrap(),
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-matrix product and sum $C =
    /// α A B + β C$ for Side is `CblasLeft` and $C = α B A + β C$ for
    /// Side is `CblasRight`, where the matrix $A$ is hermitian.  When
    /// Uplo is CblasUpper then the upper triangle and diagonal of A
    /// are used, and when Uplo is CblasLower then the lower triangle
    /// and diagonal of A are used.  The imaginary elements of the
    /// diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_zhemm")]
    pub fn hemm(
        side: Side,
        uplo: Uplo,
        alpha: &Complex<f64>,
        A: &MatrixComplexF64,
        B: &MatrixComplexF64,
        beta: &Complex<f64>,
        C: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_zhemm(
                side.into(),
                uplo.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta.unwrap(),
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes the matrix-matrix product B = \alpha
    /// op(A) B for Side is Left and B = \alpha B op(A) for Side is
    /// CblasRight.  The matrix A is triangular and op(A) = A, A^T,
    /// A^H for TransA = NoTrans, Trans, ConjTrans.  When Uplo is
    /// Upper then the upper triangle of A is used, and when Uplo is
    /// Lower then the lower triangle of A is used.  If Diag is
    /// NonUnit then the diagonal of A is used, but if Diag is Unit
    /// then the diagonal elements of the matrix A are taken as unity
    /// and are not referenced.
    #[doc(alias = "gsl_blas_ztrmm")]
    pub fn trmm(
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        alpha: &Complex<f64>,
        A: &MatrixComplexF64,
        B: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ztrmm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }

    /// This function computes the inverse-matrix matrix product B =
    /// \alpha op(inv(A))B for Side is Left and B = \alpha B
    /// op(inv(A)) for Side is Right.  The matrix A is triangular and
    /// op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and
    /// when Uplo is Lower then the lower triangle of A is used.  If
    /// Diag is NonUnit then the diagonal of A is used, but if Diag is
    /// Unit then the di agonal elements of the matrix A are taken as
    /// unity and are not referenced.
    #[doc(alias = "gsl_blas_ztrsm")]
    pub fn trsm(
        side: Side,
        uplo: Uplo,
        transA: Transpose,
        diag: Diag,
        alpha: &Complex<f64>,
        A: &MatrixComplexF64,
        B: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_ztrsm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }

    /// This function computes a rank-k update of the symmetric matrix
    /// C, C = \alpha A A^T + \beta C when Trans is NoTrans and C =
    /// \alpha A^T A + \beta C when Trans is Trans.  Since the matrix
    /// C is symmetric only its upper half or lower half need to be
    /// stored.  When Uplo is Upper then the upper triangle and
    /// diagonal of C are used, and when Uplo is Lower then the lower
    /// triangle and diagonal of C are used.
    #[doc(alias = "gsl_blas_zsyrk")]
    pub fn syrk(
        uplo: Uplo,
        trans: Transpose,
        alpha: &Complex<f64>,
        A: &MatrixComplexF64,
        beta: &Complex<f64>,
        C: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_zsyrk(
                uplo.into(),
                trans.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                beta.unwrap(),
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// These functions compute a rank-k update of the hermitian
    /// matrix C, C = \alpha A A^H + \beta C when Trans is NoTrans and
    /// C = \alpha A^H A + \beta C when Trans is ConjTrans.  Since the
    /// matrix C is hermitian only its upper half or lower half need
    /// to be stored.  When Uplo is Upper then the upper triangle and
    /// diagonal of C are used, and when Uplo is Lower then the lower
    /// triangle and diagonal of C are used.  The imaginary elements
    /// of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_zherk")]
    pub fn herk(
        uplo: Uplo,
        trans: Transpose,
        alpha: f64,
        A: &MatrixComplexF64,
        beta: f64,
        C: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_zherk(
                uplo.into(),
                trans.into(),
                alpha,
                A.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes a rank-2k update of the symmetric
    /// matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans
    /// is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when
    /// Trans is Trans.  Since the matrix C is symmetric only its
    /// upper half or lower half need to be stored.  When Uplo is
    /// Upper then the upper triangle and diagonal of C are used, and
    /// when Uplo is Lower then the lower triangle and diagonal of C
    /// are used.
    #[doc(alias = "gsl_blas_zsyr2k")]
    pub fn syr2k(
        uplo: Uplo,
        trans: Transpose,
        alpha: &Complex<f64>,
        A: &MatrixComplexF64,
        B: &MatrixComplexF64,
        beta: &Complex<f64>,
        C: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_zsyr2k(
                uplo.into(),
                trans.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta.unwrap(),
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
    /// This function computes a rank-2k update of the hermitian
    /// matrix C, C = \alpha A B^H + \alpha^* B A^H + \beta C when
    /// Trans is NoTrans and C = \alpha A^H B + \alpha^* B^H A + \beta
    /// C when Trans is ConjTrans.  Since the matrix C is hermitian
    /// only its upper half or lower half need to be stored.  When
    /// Uplo is Upper then the upper triangle and diagonal of C are
    /// used, and when Uplo is Lower then the lower triangle and
    /// diagonal of C are used.  The imaginary elements of the
    /// diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_zher2k")]
    pub fn her2k(
        uplo: Uplo,
        trans: Transpose,
        alpha: &Complex<f64>,
        A: &MatrixComplexF64,
        B: &MatrixComplexF64,
        beta: f64,
        C: &mut MatrixComplexF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_blas_zher2k(
                uplo.into(),
                trans.into(),
                alpha.unwrap(),
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_srotg() {
        let (c, s, r) = s::rotg(3., 4.).unwrap();
        assert_eq!(c, 0.6);
        assert_eq!(s, 0.8);
        assert_eq!(r, 5.);
    }

    #[test]
    fn test_drotg() {
        let (c, s, r) = d::rotg(3., 4.).unwrap();
        assert!((c - 0.6).abs() < 5e-16, "|{c} - 0.6| >= 5e-16");
        assert!((s - 0.8).abs() < 5e-16, "|{s} - 0.8| >= 5e-16");
        assert!((r - 5.).abs() < 1e-15, "|{r} - 5.| >= 1e-15");
    }
}
