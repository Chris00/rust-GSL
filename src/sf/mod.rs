//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! Special functions.
//!
//! # Usage
//!
//! The special functions are available in two calling conventions, a
//! _natural_ form which returns the numerical value of the function and
//! an _error-handling_ form which returns an error code. The two types
//! of function provide alternative ways of accessing the same
//! underlying code.

pub mod airy;
pub mod bessel;
pub mod beta;
pub mod clausen;
pub mod coulomb;
pub mod coupling_coefficients;
pub mod dawson;
pub mod debye;
pub mod dilogarithm;
pub mod elementary_operations;
pub mod elliptic;
pub mod elliptic_jacobian;
pub mod error;
pub mod exponential;
pub mod exponential_integrals;
pub mod factorials;
pub mod fermi_dirac;
pub mod gamma;
pub mod gegenbauer;
pub mod hypergeometric;
pub mod laguerre;
pub mod lambert_w;
pub mod legendre;
pub mod logarithm;
pub mod mathieu;
pub mod pochhammer_symbol;
pub mod power;
pub mod psi;
pub mod synchrotron;
pub mod transport;
pub mod trigonometric;
pub mod zeta;

/// The goal of the library is to achieve double precision accuracy
/// wherever possible.  However the cost of evaluating some special
/// functions to double precision can be significant, particularly
/// where very high order terms are required.  In these cases a mode
/// argument, of type [`Prec`] allows the accuracy of the function to
/// be reduced in order to improve performance.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Copy)]
pub enum Prec {
    /// Double-precision, a relative accuracy of approximately
    /// $2 · 10^{-16}$.
    Double,
    /// Single-precision, a relative accuracy of approximately $10^{-7}$.
    Single,
    /// Approximate values, a relative accuracy of approximately
    /// $5 · 10^{-4}$.
    Approx,
}

impl From<Prec> for sys::gsl_mode_t {
    fn from(v: Prec) -> sys::gsl_mode_t {
        match v {
            Prec::Double => sys::GSL_PREC_DOUBLE,
            Prec::Single => sys::GSL_PREC_SINGLE,
            Prec::Approx => sys::GSL_PREC_APPROX,
        }
    }
}
