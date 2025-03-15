//! Computations of the Pochhammer symbol
//! $(a)_x = \Gamma(a + x)/\Gamma(a)$.
//!
//! The gamma function is available in the [gamma module][crate::sf::gamma].

use crate::{types, Error};
use std::mem::MaybeUninit;

/// Return the Pochhammer symbol $(a)_x = \Gamma(a + x)/\Gamma(a)$.
///
/// The Pochhammer symbol is also known as the Apell symbol and
/// sometimes written as $(a,x)$.  When $a$ and $a+x$ are negative
/// integers or zero, the limiting value of the ratio is returned.
#[doc(alias = "gsl_sf_poch")]
pub fn poch(a: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_poch(a, x) }
}

/// Return the Pochhammer symbol $(a)_x = \Gamma(a + x)/\Gamma(a)$.
///
/// The Pochhammer symbol is also known as the Apell symbol and
/// sometimes written as $(a,x)$.  When $a$ and $a+x$ are negative
/// integers or zero, the limiting value of the ratio is returned.
#[doc(alias = "gsl_sf_poch_e")]
pub fn poch_e(a: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_poch_e(a, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return the logarithm of the Pochhammer symbol, $\ln((a)_x) =
/// \ln(\Gamma(a + x)/\Gamma(a))$.
#[doc(alias = "gsl_sf_lnpoch")]
pub fn lnpoch(a: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_lnpoch(a, x) }
}

/// Return the logarithm of the Pochhammer symbol, $\ln((a)_x) =
/// \ln(\Gamma(a + x)/\Gamma(a))$.
#[doc(alias = "gsl_sf_lnpoch_e")]
pub fn lnpoch_e(a: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_lnpoch_e(a, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return `(lg, sgn)` where `lg` is the logarithm of the magnitude
/// Pochhammer symbol and `sgn` its sign.
///
/// The computed parameters are `lg`$= \ln(|(a)_x|)$ with a
/// corresponding error term, and `sgn`$d= \sgn((a)_x)$ where $(a)_x =
/// \Gamma(a + x)/\Gamma(a)$.
#[doc(alias = "gsl_sf_lnpoch_sgn_e")]
pub fn lnpoch_sgn_e(a: f64, x: f64) -> Result<(types::Result, f64), Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let mut sgn = 0.;
    let ret = unsafe { sys::gsl_sf_lnpoch_sgn_e(a, x, result.as_mut_ptr(), &mut sgn) };

    Error::handle(ret, (unsafe { result.assume_init() }.into(), sgn))
}

/// Return the relative Pochhammer symbol $((a)_x - 1)/x$
/// where $(a)_x = \Gamma(a + x)/\Gamma(a)$.
#[doc(alias = "gsl_sf_pochrel")]
pub fn pochrel(a: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_pochrel(a, x) }
}

/// Return the relative Pochhammer symbol $((a)_x - 1)/x$
/// where $(a)_x = \Gamma(a + x)/\Gamma(a)$.
#[doc(alias = "gsl_sf_pochrel_e")]
pub fn pochrel_e(a: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_pochrel_e(a, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}
