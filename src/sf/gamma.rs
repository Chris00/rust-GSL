//! Compute the (incomplete) Gamma functions.
//!
//! The Gamma function is defined by the following integral,
//!
//! $$\Gamma(x) = \int_0^\infty t^{x-1} \exp(-t) dt$$
//!
//! It is related to the factorial function by $\Gamma(n)=(n-1)!$ for
//! positive integer $n$.  Further information on the Gamma function can
//! be found in Abramowitz & Stegun, Chapter 6.
//!
//! See also the modules [`factorials`][crate::sf::factorials],
//! [`pochhammer_symbol`][crate::sf::pochhammer_symbol] and
//! [`beta`][crate::sf::beta].

use crate::{types, Error};
use std::mem::MaybeUninit;

/// The maximum value of $x$ so that $Γ(x)$ does not overflow.
pub const GAMMA_XMAX: f64 = sys::GSL_SF_GAMMA_XMAX;

/// Return $Γ(x)$, subject to $x$ not being a negative integer or zero.
///
/// The function is computed using the real Lanczos method.  The
/// maximum value of $x$ such that $Γ(x)$ is not considered an
/// overflow is given by the constant [`GAMMA_XMAX`] and is 171.0.
#[doc(alias = "gsl_sf_gamma")]
pub fn gamma(x: f64) -> f64 {
    unsafe { sys::gsl_sf_gamma(x) }
}

/// Return $Γ(x)$, subject to $x$ not being a negative integer or zero.
///
/// The function is computed using the real Lanczos method.  The
/// maximum value of $x$ such that $Γ(x)$ is not considered an
/// overflow is given by the constant [`GAMMA_XMAX`] and is 171.0.
#[doc(alias = "gsl_sf_gamma_e")]
pub fn gamma_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_gamma_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\ln(Γ(x))$, subject to $x$ not being a negative integer or
/// zero.  For $x < 0$, the real part of $\ln(Γ(x))$ is returned,
/// which is equivalent to $\ln(|Γ(x)|)$.
///
/// The function is computed using the real Lanczos method.
#[doc(alias = "gsl_sf_lngamma")]
pub fn lngamma(x: f64) -> f64 {
    unsafe { sys::gsl_sf_lngamma(x) }
}

/// Return $\ln(Γ(x))$, subject to $x$ not being a negative integer or
/// zero.  For $x < 0$, the real part of $\ln(Γ(x))$ is returned,
/// which is equivalent to $\ln(|Γ(x)|)$.
///
/// The function is computed using the real Lanczos method.
#[doc(alias = "gsl_sf_lngamma_e")]
pub fn lngamma_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_lngamma_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return `(lg, sgn)` where `lg` is logarithm of the magnitude of
/// $Γ(x)$ and `sgn` is the sign of $Γ(x)$, subject to `x` not being a
/// negative integer or zero.
///
/// The function is computed using the real Lanczos method.  The value
/// of the gamma function and its error can be reconstructed using the
/// relation $Γ(x) = $`sgn` * $\exp$(`lg`), taking into account the
/// two components of `lg`.
#[doc(alias = "gsl_sf_lngamma_sgn_e")]
pub fn lngamma_sgn_e(x: f64) -> Result<(types::Result, f64), Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let mut sgn = 0.;
    let ret = unsafe { sys::gsl_sf_lngamma_sgn_e(x, result.as_mut_ptr(), &mut sgn) };

    Error::handle(ret, (unsafe { result.assume_init() }.into(), sgn))
}

/// Return the regulated Gamma function $Γ^*(x)$ for $x > 0$.
///
/// The regulated gamma function is given by,
///
/// $$Γ^*(x) = \frac{Γ(x)}{\sqrt{2\pi} x^{x-1/2} \exp(-x)}
///   = \Bigl(1 + \frac{1}{12x} + \cdots\Bigr)  \quad\text{for } x → ∞$$
///
/// and is a useful suggestion of Temme.
#[doc(alias = "gsl_sf_gammastar")]
pub fn gammastar(x: f64) -> f64 {
    unsafe { sys::gsl_sf_gammastar(x) }
}

/// Return the regulated Gamma function $Γ^*(x)$ for $x > 0$.
///
/// The regulated gamma function is given by,
///
/// $$Γ^*(x) = \frac{Γ(x)}{\sqrt{2\pi} x^{x-1/2} \exp(-x)}
///   = \Bigl(1 + \frac{1}{12x} + \cdots\Bigr)  \quad\text{for } x → ∞$$
///
/// and is a useful suggestion of Temme.
#[doc(alias = "gsl_sf_gammastar_e")]
pub fn gammastar_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_gammastar_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return the reciprocal of the gamma function, $1/Γ(x)$.
///
/// It uses the real Lanczos method.
#[doc(alias = "gsl_sf_gammainv")]
pub fn gammainv(x: f64) -> f64 {
    unsafe { sys::gsl_sf_gammainv(x) }
}

/// Return the reciprocal of the gamma function, $1/Γ(x)$.
///
/// It uses the real Lanczos method.
#[doc(alias = "gsl_sf_gammainv_e")]
pub fn gammainv_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_gammainv_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return `(lnr, arg)` being the polar representation of $\ln(Γ(z))$
/// for complex $z=$`zr` + i `zi` and $z$ not a negative integer or
/// zero.
///
/// It uses the complex Lanczos method.  The returned parameters are
/// `lnr`$= \log|Γ(z)|$ and `arg`$ = \arg(Γ(z))$ in $(-\pi,\pi]$.
/// Note that the phase part (`arg`) is not well-determined when $|z|$
/// is very large, due to inevitable roundoff in restricting to
/// $(-\pi,\pi]$.  This will result in a [`Error::Loss`] error when it
/// occurs.  The absolute value part (`lnr`), however, never suffers
/// from loss of precision.
#[doc(alias = "gsl_sf_lngamma_complex_e")]
pub fn lngamma_complex_e(zr: f64, zi: f64) -> Result<(types::Result, types::Result), Error> {
    let mut lnr = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let mut arg = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_lngamma_complex_e(zr, zi, lnr.as_mut_ptr(), arg.as_mut_ptr()) };

    // FIXME: do we want the real part to be accessible in case of
    // Error::Loss?
    Error::handle(
        ret,
        (
            unsafe { lnr.assume_init() }.into(),
            unsafe { arg.assume_init() }.into(),
        ),
    )
}

/// Return the unnormalized incomplete Gamma Function $Γ(a,x) =
/// \int_x^\infty t^{a-1} \exp(-t) dt$ for a real and $x ≥ 0$.
#[doc(alias = "gsl_sf_gamma_inc")]
pub fn gamma_inc(a: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_gamma_inc(a, x) }
}

/// Return the unnormalized incomplete Gamma Function $Γ(a,x) =
/// \int_x^\infty t^{a-1} \exp(-t) dt$ for a real and $x ≥ 0$.
#[doc(alias = "gsl_sf_gamma_inc_e")]
pub fn gamma_inc_e(a: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_gamma_inc_e(a, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return the normalized incomplete Gamma Function $Q(a,x) = 1/Γ(a)
/// \int_x^\infty t^{a-1} \exp(-t) dt$ for $a > 0$, $x ≥ 0$.
#[doc(alias = "gsl_sf_gamma_inc_Q")]
pub fn gamma_inc_Q(a: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_gamma_inc_Q(a, x) }
}

/// Return the normalized incomplete Gamma Function $Q(a,x) = 1/Γ(a)
/// \int_x^\infty t^{a-1} \exp(-t) dt$ for $a > 0$, $x ≥ 0$.
#[doc(alias = "gsl_sf_gamma_inc_Q_e")]
pub fn gamma_inc_Q_e(a: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_gamma_inc_Q_e(a, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return the complementary normalized incomplete Gamma Function
/// $P(a,x) = 1 - Q(a,x) = 1/Γ(a) \int_0^x t^{a-1} \exp(-t) dt$ for
/// $a > 0$, $x ≥ 0$.
///
/// Note that Abramowitz & Stegun call $P(a,x)$ the incomplete gamma
/// function (section 6.5).
#[doc(alias = "gsl_sf_gamma_inc_P")]
pub fn gamma_inc_P(a: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_gamma_inc_P(a, x) }
}

/// Return the complementary normalized incomplete Gamma Function
/// $P(a,x) = 1 - Q(a,x) = 1/Γ(a) \int_0^x t^{a-1} \exp(-t) dt$ for
/// $a > 0$, $x ≥ 0$.
///
/// Note that Abramowitz & Stegun call $P(a,x)$ the incomplete gamma
/// function (section 6.5).
#[doc(alias = "gsl_sf_gamma_inc_P_e")]
pub fn gamma_inc_P_e(a: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_gamma_inc_P_e(a, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}
