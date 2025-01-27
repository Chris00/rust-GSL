//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! Multiroot test algorithms, See `rgsl::types::multiroot` for solvers.

use crate::Error;
use crate::{ffi::FFI, VectorF64};

#[doc(alias = "gsl_multiroot_test_delta")]
pub fn test_delta(dx: &VectorF64, x: &VectorF64, epsabs: f64, epsrel: f64) -> Result<(), Error> {
    Error::handle(
        unsafe {
            sys::gsl_multiroot_test_delta(dx.unwrap_shared(), x.unwrap_shared(), epsabs, epsrel)
        },
        (),
    )
}

#[doc(alias = "gsl_multiroot_test_residual")]
pub fn test_residual(f: &VectorF64, epsabs: f64) -> Result<(), Error> {
    Error::handle(
        unsafe { sys::gsl_multiroot_test_residual(f.unwrap_shared(), epsabs) },
        (),
    )
}
