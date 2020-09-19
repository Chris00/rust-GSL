//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The following functions are equivalent to the function gsl_pow_int (see [1Small integer powers1](Pow.html)) with an error estimate.

use enums;
use std::mem::MaybeUninit;

/// This routine computes the power x^n for integer n. The power is computed using the minimum number of multiplications.
/// For example, x^8 is computed as ((x^2)^2)^2, requiring only 3 multiplications.
/// For reasons of efficiency, these functions do not check for overflow or underflow conditions.
///
/// ```rust
/// use rgsl::pow_int;
/// /* compute 3.0**12 */
/// println!("{}", pow_int(3., 12));
/// ```
pub fn pow_int(x: f64, n: i32) -> f64 {
    unsafe { sys::gsl_sf_pow_int(x, n) }
}

/// This routine computes the power x^n for integer n. The power is computed using the minimum number of multiplications.
/// For example, x^8 is computed as ((x^2)^2)^2, requiring only 3 multiplications.
/// For reasons of efficiency, these functions do not check for overflow or underflow conditions.
///
/// ```rust
/// use rgsl::pow_int_e;
/// /* compute 3.0**12 */
/// println!("{:?}", pow_int_e(3., 12));
/// ```
pub fn pow_int_e(x: f64, n: i32) -> Result<::types::Result, enums::Value> {
    let mut result = unsafe { MaybeUninit::<sys::gsl_sf_result>::uninit() };
    let ret = unsafe { sys::gsl_sf_pow_int_e(x, n, result.as_mut_ptr()) };

    let ret = enums::Value::from(ret);
    if ret.is_success() {
        Ok(unsafe { result.assume_init() }.into())
    } else {
        Err(ret)
    }
}
