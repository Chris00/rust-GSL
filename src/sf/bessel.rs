//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//
//! # Bessel functions
//!
//! Bessel functions are solutions to the differential equation:
//!
//! $$z^2 \frac{d^2 w}{dz^2} + z \frac{dw}{dz} + (z^2 - \nu^2) w = 0.$$

use crate::{sf::Prec, types, Error};
use std::mem::MaybeUninit;

/// Return $I_0(x)$, where $I_0$ is the regular modified cylindrical
/// Bessel function of zeroth order.
#[doc(alias = "gsl_sf_bessel_I0")]
pub fn I0(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_I0(x) }
}

/// Return $I_0(x)$, where $I_0$ is the regular modified cylindrical
/// Bessel function of zeroth order.
#[doc(alias = "gsl_sf_bessel_I0_e")]
pub fn I0_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_I0_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $I_1(x)$, where $I_1$ is the regular modified cylindrical
/// Bessel function of first order.
#[doc(alias = "gsl_sf_bessel_I1")]
pub fn I1(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_I1(x) }
}

/// Return $I_1(x)$, where $I_1$ is the regular modified cylindrical
/// Bessel function of first order.
#[doc(alias = "gsl_sf_bessel_I1_e")]
pub fn I1_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_I1_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $I_n(x)$, where $I_n$ is the regular modified cylindrical
/// Bessel function of order `n`.
#[doc(alias = "gsl_sf_bessel_In")]
pub fn In(n: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_In(n, x) }
}

/// Return $I_n(x)$, where $I_n$ is the regular modified cylindrical
/// Bessel function of order `n`.
#[doc(alias = "gsl_sf_bessel_In_e")]
pub fn In_e(n: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_In_e(n, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $I_n(x)$, where $I_n$ are the regular modified cylindrical
/// Bessel functions, for $n$ from `nmin` to `nmax` inclusive, storing
/// the results in the array `result_array`.
///
/// The start of the range `nmin` must be positive or zero.  The
/// values are computed using recurrence relations for efficiency, and
/// therefore may differ slightly from the exact values.
#[doc(alias = "gsl_sf_bessel_In_array")]
pub fn In_array(nmin: u32, nmax: u32, x: f64, result_array: &mut [f64]) -> Result<(), Error> {
    assert!(nmax - nmin < result_array.len() as _);
    let ret =
        unsafe { sys::gsl_sf_bessel_In_array(nmin as _, nmax as _, x, result_array.as_mut_ptr()) };
    Error::handle(ret, ())
}

/// Return $\exp(-|x|) I_0(x)$, the scaled regular modified
/// cylindrical Bessel function of zeroth order.
#[doc(alias = "gsl_sf_bessel_I0_scaled")]
pub fn I0_scaled(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_I0_scaled(x) }
}

/// Return $\exp(-|x|) I_0(x)$, the scaled regular modified
/// cylindrical Bessel function of zeroth order.
#[doc(alias = "gsl_sf_bessel_I0_scaled_e")]
pub fn I0_scaled_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_I0_scaled_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\exp(-|x|) I_1(x)$, the scaled regular modified
/// cylindrical Bessel function of first order.
#[doc(alias = "gsl_sf_bessel_I1_scaled")]
pub fn I1_scaled(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_I1_scaled(x) }
}

/// Return $\exp(-|x|) I_1(x)$, the scaled regular modified
/// cylindrical Bessel function of first order.
#[doc(alias = "gsl_sf_bessel_I1_scaled_e")]
pub fn I1_scaled_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_I1_scaled_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\exp(-|x|) I_n(x)$, the scaled regular modified
/// cylindrical Bessel function of order `n`.
#[doc(alias = "gsl_sf_bessel_In_scaled")]
pub fn In_scaled(n: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_In_scaled(n, x) }
}

/// Return $\exp(-|x|) I_n(x)$, the scaled regular modified
/// cylindrical Bessel function of order `n`.
#[doc(alias = "gsl_sf_bessel_In_scaled_e")]
pub fn In_scaled_e(n: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_In_scaled_e(n, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return values of the scaled regular cylindrical Bessel functions
/// $\exp(-|x|) I_n(x)$ for $n$ from `nmin` to `nmax` inclusive,
/// storing the results in the array `result_array`.
///
/// The start of the range `nmin` must be positive or zero.  The
/// values are computed using recurrence relations for efficiency, and
/// therefore may differ slightly from the exact values.
#[doc(alias = "gsl_sf_bessel_In_scaled_array")]
pub fn In_scaled_array(
    nmin: u32,
    nmax: u32,
    x: f64,
    result_array: &mut [f64],
) -> Result<(), Error> {
    assert!(nmax - nmin < result_array.len() as _);
    let ret = unsafe {
        sys::gsl_sf_bessel_In_scaled_array(nmin as _, nmax as _, x, result_array.as_mut_ptr())
    };
    Error::handle(ret, ())
}

/// Return $\exp(-|x|) i_0(x)$, the scaled regular modified spherical
/// Bessel function of zeroth order.
#[doc(alias = "gsl_sf_bessel_i0_scaled")]
pub fn i0_scaled(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_i0_scaled(x) }
}

/// Return $\exp(-|x|) i_0(x)$, the scaled regular modified spherical
/// Bessel function of zeroth order.
#[doc(alias = "gsl_sf_bessel_i0_scaled_e")]
pub fn i0_scaled_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_i0_scaled_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\exp(-|x|) i_1(x)$, the scaled regular modified spherical
/// Bessel function of first order.
#[doc(alias = "gsl_sf_bessel_i1_scaled")]
pub fn i1_scaled(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_i1_scaled(x) }
}

/// Return $\exp(-|x|) i_1(x)$, the scaled regular modified spherical
/// Bessel function of first order.
#[doc(alias = "gsl_sf_bessel_i1_scaled_e")]
pub fn i1_scaled_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_i1_scaled_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\exp(-|x|) i_2(x)$, the scaled regular modified spherical
/// Bessel function of second order.
#[doc(alias = "gsl_sf_bessel_i2_scaled")]
pub fn i2_scaled(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_i2_scaled(x) }
}

/// Return $\exp(-|x|) i_2(x)$, the scaled regular modified spherical
/// Bessel function of second order.
#[doc(alias = "gsl_sf_bessel_i2_scaled_e")]
pub fn i2_scaled_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_i2_scaled_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\exp(-|x|) i_l(x)$, the scaled regular modified spherical
/// Bessel function of order `l`.
#[doc(alias = "gsl_sf_bessel_il_scaled")]
pub fn il_scaled(l: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_il_scaled(l, x) }
}

/// Return $\exp(-|x|) i_l(x)$, the scaled regular modified spherical
/// Bessel function of order `l`.
#[doc(alias = "gsl_sf_bessel_il_scaled_e")]
pub fn il_scaled_e(l: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_il_scaled_e(l, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Store the values of the scaled regular modified cylindrical Bessel
/// functions $\exp(-|x|) i_ℓ(x)$ for $ℓ$ from $0$ to `lmax`
/// inclusive (with `lmax` ≥ 0) in the array `result_array`.
///
/// The values are computed using recurrence relations for efficiency,
/// and therefore may differ slightly from the exact values.
#[doc(alias = "gsl_sf_bessel_il_scaled_array")]
pub fn il_scaled_array(lmax: u32, x: f64, result_array: &mut [f64]) -> Result<(), Error> {
    assert!(lmax < result_array.len() as _);
    let ret =
        unsafe { sys::gsl_sf_bessel_il_scaled_array(lmax as _, x, result_array.as_mut_ptr()) };
    Error::handle(ret, ())
}

/// Return $I_\nu(x)$, the regular modified Bessel function of
/// fractional order $\nu$, for $x > 0$, $\nu > 0$.
#[doc(alias = "gsl_sf_bessel_Inu")]
pub fn Inu(nu: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_Inu(nu, x) }
}

/// Return $I_\nu(x)$, the regular modified Bessel function of
/// fractional order $\nu$, for $x > 0$, $\nu > 0$.
#[doc(alias = "gsl_sf_bessel_Inu_e")]
pub fn Inu_e(nu: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_Inu_e(nu, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\exp(-|x|)I_\nu(x)$, the scaled regular modified Bessel
/// function of fractional order $\nu$, for $x>0$, $\nu>0$.
#[doc(alias = "gsl_sf_bessel_Inu_scaled")]
pub fn Inu_scaled(nu: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_Inu_scaled(nu, x) }
}

/// Return $\exp(-|x|)I_\nu(x)$, the scaled regular modified Bessel
/// function of fractional order $\nu$, for $x>0$, $\nu>0$.
#[doc(alias = "gsl_sf_bessel_Inu_scaled_e")]
pub fn Inu_scaled_e(nu: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_Inu_scaled_e(nu, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $J_0(x)$, the regular cylindrical Bessel function of 0th
/// order.
#[doc(alias = "gsl_sf_bessel_J0")]
pub fn J0(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_J0(x) }
}

/// Return $J_0(x)$, the regular cylindrical Bessel function of 0th
/// order.
///
/// # Example
///
/// ```
/// use rgsl::sf::bessel::J0_e;
/// let x = 5.;
/// let exact = -0.17759677131433830434739701;
/// let y = J0_e(x)?;
/// println!("J₀(5) = {:.18}", y.val);
/// println!("       ± {:.18}", y.err);
/// println!("exact:  {:.18}", exact);
/// # Ok::<(), rgsl::Error>(())
/// ```
#[doc(alias = "gsl_sf_bessel_J0_e")]
pub fn J0_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_J0_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $J_1(x)$, the regular cylindrical Bessel function of first order.
#[doc(alias = "gsl_sf_bessel_J1")]
pub fn J1(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_J1(x) }
}

/// Return $J_1(x)$, the regular cylindrical Bessel function of first order.
#[doc(alias = "gsl_sf_bessel_J1_e")]
pub fn J1_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_J1_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $J_n(x)$, the regular cylindrical Bessel function of order `n`.
#[doc(alias = "gsl_sf_bessel_Jn")]
pub fn Jn(n: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_Jn(n, x) }
}

/// Return $J_n(x)$, the regular cylindrical Bessel function of order `n`.
#[doc(alias = "gsl_sf_bessel_Jn_e")]
pub fn Jn_e(n: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_Jn_e(n, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Store $J_n(x)$, the values of the regular cylindrical Bessel
/// functions, for $n$ from `nmin` to `nmax` inclusive, in the array
/// `result_array`.
///
/// The values are computed using recurrence relations for efficiency,
/// and therefore may differ slightly from the exact values.
#[doc(alias = "gsl_sf_bessel_Jn_array")]
pub fn Jn_array(nmin: u32, nmax: u32, x: f64, result_array: &mut [f64]) -> Result<(), Error> {
    assert!(nmax - nmin < result_array.len() as _);
    let ret =
        unsafe { sys::gsl_sf_bessel_Jn_array(nmin as _, nmax as _, x, result_array.as_mut_ptr()) };
    Error::handle(ret, ())
}

/// Return $j_0(x) = \sin(x)/x$, the regular spherical Bessel function
/// of zeroth order.
#[doc(alias = "gsl_sf_bessel_j0")]
pub fn j0(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_j0(x) }
}

/// Return $j_0(x) = \sin(x)/x$, the regular spherical Bessel function
/// of zeroth order.
#[doc(alias = "gsl_sf_bessel_j0_e")]
pub fn j0_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_j0_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $j_1(x) = (\sin(x)/x - \cos(x))/x$, the regular spherical
/// Bessel function of first order.
#[doc(alias = "gsl_sf_bessel_j1")]
pub fn j1(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_j1(x) }
}

/// Return $j_1(x) = (\sin(x)/x - \cos(x))/x$, the regular spherical
/// Bessel function of first order.
#[doc(alias = "gsl_sf_bessel_j1_e")]
pub fn j1_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_j1_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $j_2(x) = ((3/x^2 - 1)\sin(x) - 3\cos(x)/x)/x$, the regular
/// spherical Bessel function of second order.
#[doc(alias = "gsl_sf_bessel_j2")]
pub fn j2(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_j2(x) }
}

/// Return $j_2(x) = ((3/x^2 - 1)\sin(x) - 3\cos(x)/x)/x$, the regular
/// spherical Bessel function of second order.
#[doc(alias = "gsl_sf_bessel_j2_e")]
pub fn j2_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_j2_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $j_l(x)$, the regular spherical Bessel function of
/// order `l`, for `l` ≥ 0 and `x` ≥ 0.
#[doc(alias = "gsl_sf_bessel_jl")]
pub fn jl(l: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_jl(l, x) }
}

/// Return $j_l(x)$, the regular spherical Bessel function of
/// order `l`, for `l` ≥ 0 and `x` ≥ 0.
#[doc(alias = "gsl_sf_bessel_jl_e")]
pub fn jl_e(l: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_jl_e(l, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Store $j_ℓ(x)$, the values of the regular spherical Bessel
/// functions, for $ℓ$ from 0 to `lmax` inclusive for `lmax` ≥ 0 and
/// `x` ≥ 0, in the array `result_array`.
///
/// The values are computed using recurrence relations for efficiency,
/// and therefore may differ slightly from the exact values.
#[doc(alias = "gsl_sf_bessel_jl_array")]
pub fn jl_array(lmax: u32, x: f64, result_array: &mut [f64]) -> Result<(), Error> {
    assert!(lmax < result_array.len() as _);
    let ret = unsafe { sys::gsl_sf_bessel_jl_array(lmax as _, x, result_array.as_mut_ptr()) };
    Error::handle(ret, ())
}

/// Uses Steed’s method to compute the values of the regular spherical
/// Bessel functions $j_ℓ(x)$ for $ℓ$ from 0 to `lmax` inclusive for
/// `lmax` ≥ 0 and `x` >= 0, storing the results in the array
/// `result_array`.
///
/// The Steed/Barnett algorithm is described in Comp. Phys. Comm. 21,
/// 297 (1981).  Steed’s method is more stable than the recurrence
/// used in the other functions but is also slower.
#[doc(alias = "gsl_sf_bessel_jl_steed_array")]
pub fn jl_steed_array(lmax: u32, x: f64, result_array: &mut [f64]) -> Result<(), Error> {
    assert!(lmax < result_array.len() as _);
    let ret = unsafe { sys::gsl_sf_bessel_jl_steed_array(lmax as _, x, result_array.as_mut_ptr()) };
    Error::handle(ret, ())
}

/// Return $J_\nu(x)$, the regular cylindrical Bessel function of
/// fractional order $\nu$.
#[doc(alias = "gsl_sf_bessel_Jnu")]
pub fn Jnu(nu: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_Jnu(nu, x) }
}

/// Return $J_\nu(x)$, the regular cylindrical Bessel function of
/// fractional order $\nu$.
#[doc(alias = "gsl_sf_bessel_Jnu_e")]
pub fn Jnu_e(nu: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_Jnu_e(nu, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Evaluate $J_\nu(x)$, the regular cylindrical Bessel function of
/// fractional order $\nu$, at a series of $x$ values.
///
/// The array `v` contains the $x$ values.  They are assumed to be
/// strictly ordered and positive.  The array is over-written with the
/// values of $J_\nu(x_i)$.
#[doc(alias = "gsl_sf_bessel_sequence_Jnu_e")]
pub fn sequence_Jnu(nu: f64, mode: Prec, v: &mut [f64]) -> Result<(), Error> {
    let ret =
        unsafe { sys::gsl_sf_bessel_sequence_Jnu_e(nu, mode.into(), v.len() as _, v.as_mut_ptr()) };
    Error::handle(ret, ())
}

/// Return $K_0(x)$, the irregular modified cylindrical Bessel
/// function of zeroth order, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_K0")]
pub fn K0(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_K0(x) }
}

/// Return $K_0(x)$, the irregular modified cylindrical Bessel
/// function of zeroth order, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_K0_e")]
pub fn K0_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_K0_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $K_1(x)$, the irregular modified cylindrical Bessel
/// function of first order, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_K1")]
pub fn K1(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_K1(x) }
}

/// Return $K_1(x)$, the irregular modified cylindrical Bessel
/// function of first order, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_K1_e")]
pub fn K1_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_K1_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $K_n(x)$, the irregular modified cylindrical Bessel
/// function of order `n`, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_Kn")]
pub fn Kn(n: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_Kn(n, x) }
}

/// Return $K_n(x)$, the irregular modified cylindrical Bessel
/// function of order `n`, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_Kn_e")]
pub fn Kn_e(n: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_Kn_e(n, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Store $K_n(x)$, the values of the irregular modified cylindrical
/// Bessel functions, for $n$ from `nmin` to `nmax` inclusive, in the
/// array `result_array`.
///
/// The start of the range `nmin` must be positive or zero.  The domain
/// of the function is $x > 0$.  The values are computed using recurrence
/// relations for efficiency, and therefore may differ slightly from
/// the exact values.
#[doc(alias = "gsl_sf_bessel_Kn_array")]
pub fn Kn_array(nmin: u32, nmax: u32, x: f64, result_array: &mut [f64]) -> Result<(), Error> {
    assert!(nmax - nmin < result_array.len() as _);
    let ret =
        unsafe { sys::gsl_sf_bessel_Kn_array(nmin as _, nmax as _, x, result_array.as_mut_ptr()) };
    Error::handle(ret, ())
}

/// Return $\exp(x) K_0(x)$, the scaled irregular modified cylindrical
/// Bessel function of zeroth order for `x` > 0.
#[doc(alias = "gsl_sf_bessel_K0_scaled")]
pub fn K0_scaled(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_K0_scaled(x) }
}

/// Return $\exp(x) K_0(x)$, the scaled irregular modified cylindrical
/// Bessel function of zeroth order for `x` > 0.
#[doc(alias = "gsl_sf_bessel_K0_scaled_e")]
pub fn K0_scaled_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_K0_scaled_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\exp(x) K_1(x)$, the scaled irregular modified cylindrical
/// Bessel function of first order for `x` > 0.
#[doc(alias = "gsl_sf_bessel_K1_scaled")]
pub fn K1_scaled(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_K1_scaled(x) }
}

/// Return $\exp(x) K_1(x)$, the scaled irregular modified cylindrical
/// Bessel function of first order for `x` > 0.
#[doc(alias = "gsl_sf_bessel_K1_scaled_e")]
pub fn K1_scaled_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_K1_scaled_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\exp(x) K_n(x)$, the scaled irregular modified cylindrical
/// Bessel function of order `n`, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_Kn_scaled")]
pub fn Kn_scaled(n: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_Kn_scaled(n, x) }
}

/// Return $\exp(x) K_n(x)$, the scaled irregular modified cylindrical
/// Bessel function of order `n`, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_Kn_scaled_e")]
pub fn Kn_scaled_e(n: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_Kn_scaled_e(n, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Store $\exp(x) K_n(x)$, the values of the scaled irregular
/// cylindrical Bessel functions, for $n$ from `nmin` to `nmax`
/// inclusive, in the array `result_array`.
///
/// The start of the range `nmin` must be positive or zero.  The
/// domain of the function is $x > 0$.  The values are computed using
/// recurrence relations for efficiency, and therefore may differ
/// slightly from the exact values.
#[doc(alias = "gsl_sf_bessel_Kn_scaled_array")]
pub fn Kn_scaled_array(
    nmin: u32,
    nmax: u32,
    x: f64,
    result_array: &mut [f64],
) -> Result<(), Error> {
    assert!(nmax - nmin < result_array.len() as _);
    let ret = unsafe {
        sys::gsl_sf_bessel_Kn_scaled_array(nmin as _, nmax as _, x, result_array.as_mut_ptr())
    };
    Error::handle(ret, ())
}

/// Return $\exp(x) k_0(x)$, the scaled irregular modified spherical
/// Bessel function of zeroth order, for `x` > 0.
///
/// The irregular modified spherical Bessel functions $k_ℓ(x)$ are
/// related to the irregular modified Bessel functions of fractional
/// order, $k_ℓ(x) = \sqrt{\pi/(2x)} K_{ℓ+1/2}(x)$.
#[doc(alias = "gsl_sf_bessel_k0_scaled")]
pub fn k0_scaled(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_k0_scaled(x) }
}

/// Return $\exp(x) k_0(x)$, the scaled irregular modified spherical
/// Bessel function of zeroth order, for `x` > 0.
///
/// The irregular modified spherical Bessel functions $k_ℓ(x)$ are
/// related to the irregular modified Bessel functions of fractional
/// order, $k_ℓ(x) = \sqrt{\pi/(2x)} K_{ℓ+1/2}(x)$.
#[doc(alias = "gsl_sf_bessel_k0_scaled_e")]
pub fn k0_scaled_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_k0_scaled_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\exp(x) k_1(x),$, the scaled irregular modified spherical
/// Bessel function of first order, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_k1_scaled")]
pub fn k1_scaled(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_k1_scaled(x) }
}

/// Return $\exp(x) k_1(x),$, the scaled irregular modified spherical
/// Bessel function of first order, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_k1_scaled_e")]
pub fn k1_scaled_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_k1_scaled_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\exp(x) k_2(x)$, the scaled irregular modified spherical
/// Bessel function of second order, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_k2_scaled")]
pub fn k2_scaled(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_k2_scaled(x) }
}

/// Return $\exp(x) k_2(x)$, the scaled irregular modified spherical
/// Bessel function of second order, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_k2_scaled_e")]
pub fn k2_scaled_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_k2_scaled_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\exp(x) k_l(x)$, the scaled irregular modified spherical
/// Bessel function of order `l`, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_kl_scaled")]
pub fn kl_scaled(l: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_kl_scaled(l, x) }
}

/// Return $\exp(x) k_l(x)$, the scaled irregular modified spherical
/// Bessel function of order `l`, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_kl_scaled_e")]
pub fn kl_scaled_e(l: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_kl_scaled_e(l, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Store $\exp(x) k_ℓ(x)$, the values of the scaled irregular
/// modified spherical Bessel functions, for $ℓ$ from 0 to `lmax`
/// inclusive (for `lmax` ≥ 0 and `x` > 0), in the array `result_array`.
///
/// The values are computed using recurrence relations for efficiency,
/// and therefore may differ slightly from the exact values.
#[doc(alias = "gsl_sf_bessel_kl_scaled_array")]
pub fn kl_scaled_array(lmax: u32, x: f64, result_array: &mut [f64]) -> Result<(), Error> {
    assert!(lmax < result_array.len() as _);
    let ret =
        unsafe { sys::gsl_sf_bessel_kl_scaled_array(lmax as _, x, result_array.as_mut_ptr()) };
    Error::handle(ret, ())
}

/// Return $K_\nu(x)$, the irregular modified Bessel function of
/// fractional order $\nu$, for $x>0$, $\nu>0$.
#[doc(alias = "gsl_sf_bessel_Knu")]
pub fn Knu(nu: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_Knu(nu, x) }
}

/// Return $K_\nu(x)$, the irregular modified Bessel function of
/// fractional order $\nu$, for $x>0$, $\nu>0$.
#[doc(alias = "gsl_sf_bessel_Knu_e")]
pub fn Knu_e(nu: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_Knu_e(nu, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\ln(K_\nu(x))$, the logarithm of the irregular modified
/// Bessel function of fractional order $\nu$, for $x>0$, $\nu>0$.
#[doc(alias = "gsl_sf_bessel_lnKnu")]
pub fn lnKnu(nu: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_lnKnu(nu, x) }
}

/// Return $\ln(K_\nu(x))$, the logarithm of the irregular modified
/// Bessel function of fractional order $\nu$, for $x>0$, $\nu>0$.
#[doc(alias = "gsl_sf_bessel_lnKnu_e")]
pub fn lnKnu_e(nu: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_lnKnu_e(nu, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $\exp(+|x|) K_\nu(x)$, the scaled irregular modified Bessel
/// function of fractional order $\nu$, for $x>0$, $\nu>0$.
#[doc(alias = "gsl_sf_bessel_Knu_scaled")]
pub fn Knu_scaled(nu: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_Knu_scaled(nu, x) }
}

/// Return $\exp(+|x|) K_\nu(x)$, the scaled irregular modified Bessel
/// function of fractional order $\nu$, for $x>0$, $\nu>0$.
#[doc(alias = "gsl_sf_bessel_Knu_scaled_e")]
pub fn Knu_scaled_e(nu: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_Knu_scaled_e(nu, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $Y_0(x)$, the irregular cylindrical Bessel function of
/// zeroth order, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_Y0")]
pub fn Y0(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_Y0(x) }
}

/// Return $Y_0(x)$, the irregular cylindrical Bessel function of
/// zeroth order, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_Y0_e")]
pub fn Y0_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_Y0_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $Y_1(x)$, the irregular cylindrical Bessel function of
/// first order, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_Y1")]
pub fn Y1(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_Y1(x) }
}

/// Return $Y_1(x)$, the irregular cylindrical Bessel function of
/// first order, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_Y1_e")]
pub fn Y1_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_Y1_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $Y_n(x)$, the irregular cylindrical Bessel function of
/// order `n`, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_Yn")]
pub fn Yn(n: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_Yn(n, x) }
}

/// Return $Y_n(x)$, the irregular cylindrical Bessel function of
/// order `n`, for `x` > 0.
#[doc(alias = "gsl_sf_bessel_Yn_e")]
pub fn Yn_e(n: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_Yn_e(n, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Store $Y_n(x)$, the values of the irregular cylindrical Bessel
/// functions, for $n$ from `nmin` to `nmax` inclusive, in the array
/// `result_array`.
///
/// The domain of the function is $x>0$.  The values are computed
/// using recurrence relations for efficiency, and therefore may
/// differ slightly from the exact values.
#[doc(alias = "gsl_sf_bessel_Yn_array")]
pub fn Yn_array(nmin: u32, nmax: u32, x: f64, result_array: &mut [f64]) -> Result<(), Error> {
    assert!(nmax - nmin < result_array.len() as _);
    let ret =
        unsafe { sys::gsl_sf_bessel_Yn_array(nmin as _, nmax as _, x, result_array.as_mut_ptr()) };
    Error::handle(ret, ())
}

/// Return $y_0(x) = -\cos(x)/x$, the irregular spherical Bessel
/// function of zeroth order.
#[doc(alias = "gsl_sf_bessel_y0")]
pub fn y0(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_y0(x) }
}

/// Return $y_0(x) = -\cos(x)/x$, the irregular spherical Bessel
/// function of zeroth order.
#[doc(alias = "gsl_sf_bessel_y0_e")]
pub fn y0_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_y0_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $y_1(x) = -(\cos(x)/x + \sin(x))/x$, the irregular
/// spherical Bessel function of first order.
#[doc(alias = "gsl_sf_bessel_y1")]
pub fn y1(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_y1(x) }
}

/// Return $y_1(x) = -(\cos(x)/x + \sin(x))/x$, the irregular
/// spherical Bessel function of first order.
#[doc(alias = "gsl_sf_bessel_y1_e")]
pub fn y1_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_y1_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $y_2(x) = (-3/x^3 + 1/x)\cos(x) - (3/x^2)\sin(x)$, the
/// irregular spherical Bessel function of second order.
#[doc(alias = "gsl_sf_bessel_y2")]
pub fn y2(x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_y2(x) }
}

/// Return $y_2(x) = (-3/x^3 + 1/x)\cos(x) - (3/x^2)\sin(x)$, the
/// irregular spherical Bessel function of second order.
#[doc(alias = "gsl_sf_bessel_y2_e")]
pub fn y2_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_y2_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return $y_l(x)$, the irregular spherical Bessel function of order
/// `l`, for `l` ≥ 0.
#[doc(alias = "gsl_sf_bessel_yl")]
pub fn yl(l: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_yl(l, x) }
}

/// Return $y_l(x)$, the irregular spherical Bessel function of order
/// `l`, for `l` ≥ 0.
#[doc(alias = "gsl_sf_bessel_yl_e")]
pub fn yl_e(l: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_yl_e(l, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Store $y_ℓ(x)$, the values of the irregular spherical Bessel
/// functions, for $ℓ$ from 0 to `lmax` inclusive for `lmax` ≥ 0, in
/// the array `result_array`.
///
/// The values are computed using recurrence relations for efficiency,
/// and therefore may differ slightly from the exact values.
#[doc(alias = "gsl_sf_bessel_yl_array")]
pub fn yl_array(lmax: u32, x: f64, result_array: &mut [f64]) -> Result<(), Error> {
    assert!(lmax < result_array.len() as _);
    let ret = unsafe { sys::gsl_sf_bessel_yl_array(lmax as _, x, result_array.as_mut_ptr()) };
    Error::handle(ret, ())
}

/// Return $Y_\nu(x)$, the irregular cylindrical Bessel function of
/// fractional order $\nu$.
#[doc(alias = "gsl_sf_bessel_Ynu")]
pub fn Ynu(nu: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_bessel_Ynu(nu, x) }
}

/// Return $Y_\nu(x)$, the irregular cylindrical Bessel function of
/// fractional order $\nu$.
#[doc(alias = "gsl_sf_bessel_Ynu_e")]
pub fn Ynu_e(nu: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_Ynu_e(nu, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return the location of the `s`-th positive zero of the Bessel
/// function $J_0(x)$.
#[doc(alias = "gsl_sf_bessel_zero_J0")]
pub fn zero_J0(s: u32) -> f64 {
    unsafe { sys::gsl_sf_bessel_zero_J0(s) }
}

/// Return the location of the `s`-th positive zero of the Bessel
/// function $J_0(x)$.
#[doc(alias = "gsl_sf_bessel_zero_J0_e")]
pub fn zero_J0_e(s: u32) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_zero_J0_e(s, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return the location of the `s`-th positive zero of the Bessel
/// function $J_1(x)$.
#[doc(alias = "gsl_sf_bessel_zero_J1")]
pub fn zero_J1(s: u32) -> f64 {
    unsafe { sys::gsl_sf_bessel_zero_J1(s) }
}

/// Return the location of the `s`-th positive zero of the Bessel
/// function $J_1(x)$.
#[doc(alias = "gsl_sf_bessel_zero_J1_e")]
pub fn zero_J1_e(s: u32) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_zero_J1_e(s, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// Return the location of the `s`-th positive zero of the Bessel
/// function $J_\nu(x)$.
///
/// The current implementation does not support negative values of $\nu$.
#[doc(alias = "gsl_sf_bessel_zero_Jnu")]
pub fn zero_Jnu(nu: f64, s: u32) -> f64 {
    unsafe { sys::gsl_sf_bessel_zero_Jnu(nu, s) }
}

/// Return the location of the `s`-th positive zero of the Bessel
/// function $J_\nu(x)$.
///
/// The current implementation does not support negative values of $\nu$.
#[doc(alias = "gsl_sf_bessel_zero_Jnu_e")]
pub fn zero_Jnu_e(nu: f64, s: u32) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_bessel_zero_Jnu_e(nu, s, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}
