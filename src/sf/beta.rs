//! Beta functions.
//!
//! The Gamma function is available in the module
//! [`sf::gamma`][crate::sf::gamma].

use crate::{types, Error};
use std::mem::MaybeUninit;

/// This routine computes the Beta Function, $B(a,b) =
/// \Gamma(a)\Gamma(b)/\Gamma(a+b)$ subject to $a$ and $b$ not being
/// negative integers.
#[doc(alias = "gsl_sf_beta")]
pub fn beta(a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_sf_beta(a, b) }
}

/// This routine computes the Beta Function, $B(a,b) =
/// \Gamma(a)\Gamma(b)/\Gamma(a+b)$ subject to $a$ and $b$ not being
/// negative integers.
#[doc(alias = "gsl_sf_beta_e")]
pub fn beta_e(a: f64, b: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_beta_e(a, b, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the logarithm of the Beta Function,
/// $\log(B(a,b))$ subject to $a$ and $b$ not being negative integers.
#[doc(alias = "gsl_sf_lnbeta")]
pub fn lnbeta(a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_sf_lnbeta(a, b) }
}

/// This routine computes the logarithm of the Beta Function,
/// $\log(B(a,b))$ subject to $a$ and $b$ not being negative integers.
#[doc(alias = "gsl_sf_lnbeta_e")]
pub fn lnbeta_e(a: f64, b: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_lnbeta_e(a, b, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the normalized incomplete Beta function
/// $I_x(a,b) = B_x(a,b)/B(a,b)$.
///
/// The function $B_x$ is defined by
///
/// $$B_x(a,b) = \int_0^x t^{a-1} (1-t)^{b-1} dt$$
///
/// for $0 ≤ x ≤ 1$.
///
/// For $a > 0$, $b > 0$ the value is computed using a continued
/// fraction expansion.  For all other values it is computed using the
/// relation $I_x(a,b,x) = (1/a) x^a 2F1(a,1-b,a+1,x)/B(a,b)$.
#[doc(alias = "gsl_sf_beta_inc")]
pub fn beta_inc(a: f64, b: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_beta_inc(a, b, x) }
}

/// This routine computes the normalized incomplete Beta function
/// $I_x(a,b)=B_x(a,b)/B(a,b)$.
///
/// The function $B_x$ is defined by
///
/// $$B_x(a,b) = \int_0^x t^{a-1} (1-t)^{b-1} dt$$
///
/// for $0 ≤ x ≤ 1$.
///
/// For $a > 0$, $b > 0$ the value is computed using a continued
/// fraction expansion.  For all other values it is computed using the
/// relation $I_x(a,b,x) = (1/a) x^a 2F1(a,1-b,a+1,x)/B(a,b)$.
#[doc(alias = "gsl_sf_beta_inc_e")]
pub fn beta_inc_e(a: f64, b: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_beta_inc_e(a, b, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}
