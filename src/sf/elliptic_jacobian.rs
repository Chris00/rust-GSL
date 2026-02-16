//! The Jacobian Elliptic functions are defined in Abramowitz &
//! Stegun, Chapter 16.

use crate::Error;

/// The values of the Jacobian elliptic functions $\sn(u|m)$,
/// $\cn(u|m)$, and $\dn(u|m)$.
#[derive(Clone, Copy, Debug)]
pub struct EllJac {
    pub sn: f64,
    pub cn: f64,
    pub dn: f64,
}

/// Return the values of the Jacobian elliptic functions $\sn(u|m)$,
/// $\cn(u|m)$, $\dn(u|m)$.
///
/// It uses descending Landen transformations.
#[doc(alias = "gsl_sf_elljac_e")]
pub fn elljac_e(u: f64, m: f64) -> Result<EllJac, Error> {
    let mut sn = 0.;
    let mut cn = 0.;
    let mut dn = 0.;
    let ret = unsafe { sys::gsl_sf_elljac_e(u, m, &mut sn, &mut cn, &mut dn) };
    Error::handle(ret, EllJac { sn, cn, dn })
}

/// Return $\sn(u|m)$.
///
/// This is a convenience function.  It you need to compute other
/// Jacobian elliptic functions for the same `(u, m)`, use [`elljac_e`].
pub fn sn(u: f64, m: f64) -> Result<f64, Error> {
    elljac_e(u, m).map(|e| e.sn)
}

/// Return $\cn(u|m)$.
///
/// This is a convenience function.  It you need to compute other
/// Jacobian elliptic functions for the same `(u, m)`, use [`elljac_e`].
pub fn cn(u: f64, m: f64) -> Result<f64, Error> {
    elljac_e(u, m).map(|e| e.cn)
}

/// Return $\dn(u|m)$.
///
/// This is a convenience function.  It you need to compute other
/// Jacobian elliptic functions for the same `(u, m)`, use [`elljac_e`].
pub fn dn(u: f64, m: f64) -> Result<f64, Error> {
    elljac_e(u, m).map(|e| e.dn)
}
