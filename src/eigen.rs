//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*! Eigenvalues and eigenvectors.

## References and Further Reading

Further information on the algorithms described in this section can be found in the following book,

G. H. Golub, C. F. Van Loan, Matrix Computations (3rd Ed, 1996), Johns Hopkins University Press, ISBN 0-8018-5414-8.
Further information on the generalized eigensystems QZ algorithm can be found in this paper,

C. Moler, G. Stewart, “An Algorithm for Generalized Matrix Eigenvalue Problems”, SIAM J. Numer. Anal., Vol 10, No 2, 1973.
Eigensystem routines for very large matrices can be found in the Fortran library LAPACK. The LAPACK library is described in,

LAPACK Users’ Guide (Third Edition, 1999), Published by SIAM, ISBN 0-89871-447-8.
<http://www.netlib.org/lapack>

The LAPACK source code can be found at the website above along with an online copy of the users guide.
*/

use crate::ffi::FFI;
#[cfg(feature = "complex")]
use crate::types::{MatrixComplexF64, VectorComplexF64};
use crate::types::{MatrixF64, VectorF64};
use crate::Error;

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Copy)]
pub enum Sort {
    /// ascending order in numerical value
    Asc,
    /// descending order in numerical value
    Desc,
    /// ascending order in magnitude
    AbsAsc,
    /// descending order in magnitude
    AbsDesc,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::gsl_eigen_sort_t> for Sort {
    fn into(self) -> sys::gsl_eigen_sort_t {
        match self {
            Self::Asc => sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_VAL_ASC,
            Self::Desc => sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_VAL_DESC,
            Self::AbsAsc => sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_ABS_ASC,
            Self::AbsDesc => sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_ABS_DESC,
        }
    }
}

#[doc(hidden)]
impl From<sys::gsl_eigen_sort_t> for Sort {
    fn from(v: sys::gsl_eigen_sort_t) -> Sort {
        match v {
            sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_VAL_ASC => Self::Asc,
            sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_VAL_DESC => Self::Desc,
            sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_ABS_ASC => Self::AbsAsc,
            sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_ABS_DESC => Self::AbsDesc,
            _ => panic!("Unknown rgsl::eigen::Sort value"),
        }
    }
}

/// This function simultaneously sorts the eigenvalues stored in the
/// vector eval and the corresponding real eigenvectors stored in the
/// columns of the matrix evec into ascending or descending order
/// according to the value of the parameter `sort`.
#[doc(alias = "gsl_eigen_symmv_sort")]
pub fn symmv_sort(eval: &mut VectorF64, evec: &mut MatrixF64, sort: Sort) -> Result<(), Error> {
    let ret = unsafe {
        sys::gsl_eigen_symmv_sort(eval.unwrap_unique(), evec.unwrap_unique(), sort.into())
    };
    Error::handle(ret, ())
}

/// This function simultaneously sorts the eigenvalues stored in the
/// vector eval and the corresponding complex eigenvectors stored in
/// the columns of the matrix evec into ascending or descending order
/// according to the value of the parameter `sort`.
#[doc(alias = "gsl_eigen_hermv_sort")]
#[cfg(feature = "complex")]
pub fn hermv_sort(
    eval: &mut VectorF64,
    evec: &mut MatrixComplexF64,
    sort: Sort,
) -> Result<(), Error> {
    let ret = unsafe {
        sys::gsl_eigen_hermv_sort(eval.unwrap_unique(), evec.unwrap_unique(), sort.into())
    };
    Error::handle(ret, ())
}

/// This function simultaneously sorts the eigenvalues stored in the
/// vector eval and the corresponding complex eigenvectors stored in
/// the columns of the matrix evec into ascending or descending order
/// according to the value of the parameter `sort`.  Only
/// [`Sort::AbsAsc`] and [`Sort::AbsDesc`] are supported due to the
/// eigenvalues being complex.
#[doc(alias = "gsl_eigen_nonsymmv_sort")]
#[cfg(feature = "complex")]
pub fn nonsymmv_sort(
    eval: &mut VectorComplexF64,
    evec: &mut MatrixComplexF64,
    sort: Sort,
) -> Result<(), Error> {
    let ret = unsafe {
        sys::gsl_eigen_nonsymmv_sort(eval.unwrap_unique(), evec.unwrap_unique(), sort.into())
    };
    Error::handle(ret, ())
}

/// This function simultaneously sorts the eigenvalues stored in the
/// vector eval and the corresponding real eigenvectors stored in the
/// columns of the matrix evec into ascending or descending order
/// according to the value of the parameter `sort`.
#[doc(alias = "gsl_eigen_gensymmv_sort")]
pub fn gensymmv_sort(eval: &mut VectorF64, evec: &mut MatrixF64, sort: Sort) -> Result<(), Error> {
    let ret = unsafe {
        sys::gsl_eigen_gensymmv_sort(eval.unwrap_unique(), evec.unwrap_unique(), sort.into())
    };
    Error::handle(ret, ())
}

/// This function simultaneously sorts the eigenvalues stored in the
/// vector eval and the corresponding complex eigenvectors stored in
/// the columns of the matrix evec into ascending or descending order
/// according to the value of the parameter `sort`.
#[doc(alias = "gsl_eigen_genhermv_sort")]
#[cfg(feature = "complex")]
pub fn genhermv_sort(
    eval: &mut VectorF64,
    evec: &mut MatrixComplexF64,
    sort: Sort,
) -> Result<(), Error> {
    let ret = unsafe {
        sys::gsl_eigen_genhermv_sort(eval.unwrap_unique(), evec.unwrap_unique(), sort.into())
    };
    Error::handle(ret, ())
}

/// This function simultaneously sorts the eigenvalues stored in the
/// vectors (alpha, beta) and the corresponding complex eigenvectors
/// stored in the columns of the matrix evec into ascending or
/// descending order according to the value of the parameter `sort`.
/// Only [`Sort::AbsAsc`] and [`Sort::AbsDesc`] are supported due to
/// the eigenvalues being complex.
#[doc(alias = "gsl_eigen_genv_sort")]
#[cfg(feature = "complex")]
pub fn genv_sort(
    alpha: &mut VectorComplexF64,
    beta: &mut VectorF64,
    evec: &mut MatrixComplexF64,
    sort: Sort,
) -> Result<(), Error> {
    let ret = unsafe {
        sys::gsl_eigen_genv_sort(
            alpha.unwrap_unique(),
            beta.unwrap_unique(),
            evec.unwrap_unique(),
            sort.into(),
        )
    };
    Error::handle(ret, ())
}
