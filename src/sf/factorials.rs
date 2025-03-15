//! Although factorials can be computed from the
//! [`gamma`][crate::sf::gamma::gamma] function, using the relation
//! $n!  = \Gamma(n+1)$ for non-negative integer $n$, it is usually
//! more efficient to call the functions in this module, particularly
//! for small values of $n$, whose factorial values are maintained in
//! hardcoded tables.

use crate::{types, Error};
use std::mem::MaybeUninit;

/// This routine computes the factorial n!. The factorial is related
/// to the Gamma function by n! = \Gamma(n+1).  The maximum value of n
/// such that n! is not considered an overflow is given by the macro
/// SF_FACT_NMAX and is 170.
#[doc(alias = "gsl_sf_fact")]
pub fn fact(n: u32) -> f64 {
    unsafe { sys::gsl_sf_fact(n) }
}

/// This routine computes the factorial n!. The factorial is related
/// to the Gamma function by n! = \Gamma(n+1).  The maximum value of n
/// such that n! is not considered an overflow is given by the macro
/// SF_FACT_NMAX and is 170.
#[doc(alias = "gsl_sf_fact_e")]
pub fn fact_e(n: u32) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_fact_e(n, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the double factorial n!! = n(n-2)(n-4)
/// \dots.  The maximum value of n such that n!! is not considered an
/// overflow is given by the macro SF_DOUBLEFACT_NMAX and is 297.
#[doc(alias = "gsl_sf_doublefact")]
pub fn doublefact(n: u32) -> f64 {
    unsafe { sys::gsl_sf_doublefact(n) }
}

/// This routine computes the double factorial n!! = n(n-2)(n-4)
/// \dots.  The maximum value of n such that n!! is not considered an
/// overflow is given by the macro SF_DOUBLEFACT_NMAX and is 297.
#[doc(alias = "gsl_sf_doublefact_e")]
pub fn doublefact_e(n: u32) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_doublefact_e(n, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the logarithm of the factorial of n,
/// \log(n!).  The algorithm is faster than computing \ln(\Gamma(n+1))
/// via gsl_sf_lngamma for n < 170, but defers for larger n.
#[doc(alias = "gsl_sf_lnfact")]
pub fn lnfact(n: u32) -> f64 {
    unsafe { sys::gsl_sf_lnfact(n) }
}

/// This routine computes the logarithm of the factorial of n,
/// \log(n!).  The algorithm is faster than computing \ln(\Gamma(n+1))
/// via gsl_sf_lngamma for n < 170, but defers for larger n.
#[doc(alias = "gsl_sf_lnfact_e")]
pub fn lnfact_e(n: u32) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_lnfact_e(n, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the logarithm of the double factorial of n,
/// \log(n!!).
#[doc(alias = "gsl_sf_lndoublefact")]
pub fn lndoublefact(n: u32) -> f64 {
    unsafe { sys::gsl_sf_lndoublefact(n) }
}

/// This routine computes the logarithm of the double factorial of n,
/// \log(n!!).
#[doc(alias = "gsl_sf_lndoublefact_e")]
pub fn lndoublefact_e(n: u32) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_lndoublefact_e(n, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the combinatorial factor n choose m = n!/(m!(n-m)!)
#[doc(alias = "gsl_sf_choose")]
pub fn choose(n: u32, m: u32) -> f64 {
    unsafe { sys::gsl_sf_choose(n, m) }
}

/// This routine computes the combinatorial factor n choose m = n!/(m!(n-m)!)
#[doc(alias = "gsl_sf_choose_e")]
pub fn choose_e(n: u32, m: u32) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_choose_e(n, m, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the logarithm of n choose m. This is
/// equivalent to the sum \log(n!) - \log(m!) - \log((n-m)!).
#[doc(alias = "gsl_sf_lnchoose")]
pub fn lnchoose(n: u32, m: u32) -> f64 {
    unsafe { sys::gsl_sf_lnchoose(n, m) }
}

/// This routine computes the logarithm of n choose m. This is
/// equivalent to the sum \log(n!) - \log(m!) - \log((n-m)!).
#[doc(alias = "gsl_sf_lnchoose_e")]
pub fn lnchoose_e(n: u32, m: u32) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_lnchoose_e(n, m, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the Taylor coefficient x^n / n! for x >= 0, n >= 0.
#[doc(alias = "gsl_sf_taylorcoeff")]
pub fn taylorcoeff(n: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_taylorcoeff(n, x) }
}

/// This routine computes the Taylor coefficient x^n / n! for x >= 0, n >= 0.
#[doc(alias = "gsl_sf_taylorcoeff_e")]
pub fn taylorcoeff_e(n: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_taylorcoeff_e(n, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}
