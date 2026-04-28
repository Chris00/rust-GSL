//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::ffi::FFI;
use crate::{Error, MatrixF64, vector::VecF64};

/// Compute the covariance matrix cov = inv (J^T J) by QRP^T decomposition of J
#[doc(alias = "gsl_multifit_covar")]
pub fn covar(J: &MatrixF64, epsrel: f64, covar: &mut MatrixF64) -> Result<(), Error> {
    let ret = unsafe { sys::gsl_multifit_covar(J.unwrap_shared(), epsrel, covar.unwrap_unique()) };
    Error::handle(ret, ())
}

#[doc(alias = "gsl_multifit_test_delta")]
pub fn test_delta(dx: &VecF64, x: &VecF64, epsabs: f64, epsrel: f64) -> Result<(), Error> {
    let ret = unsafe {
        sys::gsl_multifit_test_delta(dx.unwrap_shared(), x.unwrap_shared(), epsabs, epsrel)
    };
    Error::handle(ret, ())
}

#[doc(alias = "gsl_multifit_gradient")]
pub fn gradient(J: &MatrixF64, f: &VecF64, g: &mut VecF64) -> Result<(), Error> {
    let ret = unsafe {
        sys::gsl_multifit_gradient(J.unwrap_shared(), f.unwrap_shared(), g.unwrap_unique())
    };
    Error::handle(ret, ())
}

#[doc(alias = "gsl_multifit_linear_lreg")]
pub fn linear_lreg(smin: f64, smax: f64, reg_param: &mut VecF64) -> Result<(), Error> {
    let ret = unsafe { sys::gsl_multifit_linear_lreg(smin, smax, reg_param.unwrap_unique()) };
    Error::handle(ret, ())
}

/// Returns `idx`.
#[doc(alias = "gsl_multifit_linear_lcorner")]
pub fn linear_lcorner(rho: &VecF64, eta: &VecF64) -> Result<usize, Error> {
    let mut idx = 0;
    let ret = unsafe {
        sys::gsl_multifit_linear_lcorner(rho.unwrap_shared(), eta.unwrap_shared(), &mut idx)
    };
    Error::handle(ret, idx)
}

/// Returns `(Value, idx)`.
#[doc(alias = "gsl_multifit_linear_lcorner2")]
pub fn linear_lcorner2(rho: &VecF64, eta: &VecF64) -> Result<usize, Error> {
    let mut idx = 0;
    let ret = unsafe {
        sys::gsl_multifit_linear_lcorner2(rho.unwrap_shared(), eta.unwrap_shared(), &mut idx)
    };
    Error::handle(ret, idx)
}

#[doc(alias = "gsl_multifit_linear_Lk")]
pub fn linear_Lk(p: usize, k: usize, L: &mut MatrixF64) -> Result<(), Error> {
    let ret = unsafe { sys::gsl_multifit_linear_Lk(p, k, L.unwrap_unique()) };
    Error::handle(ret, ())
}
