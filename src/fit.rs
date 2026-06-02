//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Linear Regression

This module contains routines for performing least squares fits to
experimental data using linear combinations of functions.  The data
may be weighted or unweighted, i.e. with known or unknown errors.  For
weighted data the functions compute the best fit parameters and their
associated covariance matrix. For unweighted data the covariance
matrix is estimated from the scatter of the points, giving a
variance-covariance matrix.

The functions are divided into separate versions for simple one- or
two-parameter regression and multiple-parameter fits.

# Overview

Least-squares fits are found by minimizing $χ²$ (chi-squared), the
weighted sum of squared residuals over $n$ experimental datapoints
$(xᵢ, yᵢ)$ for the model $Y(c,x)$,

$$χ² = ∑ᵢ wᵢ (yᵢ - Y(c, xᵢ))²$$

The $p$ parameters of the model are $c = \{c_0, c_1,…\}$.  The weight
factors $wᵢ$ are given by $wᵢ = 1/σᵢ²$ where $σᵢ$ is the experimental
error on the data-point $yᵢ$.  The errors are assumed to be Gaussian
and uncorrelated.  For unweighted data the chi-squared sum is computed
without any weight factors.

The fitting routines return the best-fit parameters $c$ and their $p ×
p$ covariance matrix.  The covariance matrix measures the statistical
errors on the best-fit parameters resulting from the errors on the
data, $σᵢ$, and is defined as

$$C_{ab} = ⟨δcₐ δc_b⟩

where $⟨\,⟩$ denotes an average over the Gaussian error distributions
of the underlying datapoints.

The covariance matrix is calculated by error propagation from the data
errors $σᵢ$.  The change in a fitted parameter $δcₐ$ caused by a small
change in the data $δyᵢ$ is given by

$$δcₐ = ∑ᵢ \frac{∂cₐ}{∂yᵢ} δyᵢ$$

allowing the covariance matrix to be written in terms of the errors on
the data,

$$C_{ab} = ∑_{i,j} \frac{∂cₐ}{∂yᵢ} \frac{∂c_b}{∂yⱼ} ⟨δyᵢ δyⱼ⟩$$

For uncorrelated data the fluctuations of the underlying datapoints satisfy

$$⟨δyᵢ δyⱼ⟩ = σᵢ² δ_{ij}$$

giving a corresponding parameter covariance matrix of

$$C_{ab} = ∑ᵢ \frac{1}{wᵢ} \frac{∂cₐ}{∂yᵢ} \frac{∂c_b}{∂yᵢ}$$

When computing the covariance matrix for unweighted data, i.e. data
with unknown errors, the weight factors $wᵢ$ in this sum are replaced
by the single estimate $w = 1/σ²$, where $|sigma²$ is the computed
variance of the residuals about the best-fit model, $σ² = ∑ (yᵢ -
Y(c,xᵢ))² / (n-p)$.  This is referred to as the variance-covariance
matrix.

The standard deviations of the best-fit parameters are given by the
square root of the corresponding diagonal elements of the covariance
matrix, $σ_{c_a} = √{C_{aa}}$.  The correlation coefficient of the fit
parameters $c_a$ and $c_b$ is given by $ρ_{ab} = C_{ab} / √{C_{aa}
C_{bb}}$.

*/

use crate::{
    Error,
    vector::{Vector, check_equal_len},
};

/// This function computes the best-fit linear regression coefficients
/// (c0, c1) of the model Y = c_0 + c_1 X for the dataset (`x`, `y`),
/// two vectors of the same length (possibly with strides).
///
/// The errors on `y` are assumed unknown so the variance-covariance
/// matrix for the parameters (c0, c1) is estimated from the scatter
/// of the points around the best-fit line and returned via the
/// parameters (cov00, cov01, cov11).
///
/// The sum of squares of the residuals from the best-fit line is
/// returned in sumsq. Note: the correlation coefficient of the data
/// can be computed using gsl_stats_correlation (see
/// [`Correlation`](http://www.gnu.org/software/gsl/manual/html_node/Correlation.html#Correlation)),
/// it does not depend on the fit.
///
/// Returns `(c0, c1, cov00, cov01, cov11, sumsq)`.
///
/// # Example
///
/// ```
/// use rgsl::fit;
/// let (c0, c1, _, _, _, _) = fit::linear(&[0., 1.], &[0., 1.])?;
/// assert_eq!(c0, 0.);
/// assert_eq!(c1, 1.);
/// # Ok::<(), rgsl::Error>(())
/// ```
#[doc(alias = "gsl_fit_linear")]
pub fn linear<T>(x: &T, y: &T) -> Result<(f64, f64, f64, f64, f64, f64), Error>
where
    T: Vector<f64> + ?Sized,
{
    check_equal_len(x, y)?;
    let mut c0 = 0.;
    let mut c1 = 0.;
    let mut cov00 = 0.;
    let mut cov01 = 0.;
    let mut cov11 = 0.;
    let mut sumsq = 0.;
    let ret = unsafe {
        ::sys::gsl_fit_linear(
            T::as_slice(x).as_ptr(),
            T::stride(x),
            T::as_slice(y).as_ptr(),
            T::stride(y),
            T::len(x),
            &mut c0,
            &mut c1,
            &mut cov00,
            &mut cov01,
            &mut cov11,
            &mut sumsq,
        )
    };
    Error::handle(ret, (c0, c1, cov00, cov01, cov11, sumsq))
}

/// This function computes the best-fit linear regression coefficients
/// (c0, c1) of the model Y = c_0 + c_1 X for the weighted dataset
/// (`x`, `y`), two vectors of the same length (possibly with strides).
///
/// The vector `w`, of the same length as `x` and `y`, specifies the
/// weight of each datapoint.
///
/// The weight is the reciprocal of the variance for each datapoint in y.
///
/// The covariance matrix for the parameters (c0, c1) is computed using the weights and returned via
/// the parameters (cov00, cov01, cov11).
/// The weighted sum of squares of the residuals from the best-fit line, \chi^2, is returned in chisq.
///
/// Returns `(c0, c1, cov00, cov01, cov11, chisq)`.
#[doc(alias = "gsl_fit_wlinear")]
pub fn wlinear<T: Vector<f64> + ?Sized>(
    x: &T,
    w: &T,
    y: &T,
) -> Result<(f64, f64, f64, f64, f64, f64), Error> {
    check_equal_len(x, y)?;
    check_equal_len(x, w)?;
    let mut c0 = 0.;
    let mut c1 = 0.;
    let mut cov00 = 0.;
    let mut cov01 = 0.;
    let mut cov11 = 0.;
    let mut chisq = 0.;
    let ret = unsafe {
        ::sys::gsl_fit_wlinear(
            T::as_slice(x).as_ptr(),
            T::stride(x),
            T::as_slice(w).as_ptr(),
            T::stride(w),
            T::as_slice(y).as_ptr(),
            T::stride(y),
            T::len(x),
            &mut c0,
            &mut c1,
            &mut cov00,
            &mut cov01,
            &mut cov11,
            &mut chisq,
        )
    };
    Error::handle(ret, (c0, c1, cov00, cov01, cov11, chisq))
}

/// This function uses the best-fit linear regression coefficients c0, c1 and their covariance
/// cov00, cov01, cov11 to compute the fitted function y and its standard deviation y_err for the
/// model Y = c_0 + c_1 X at the point x.
///
/// Returns `(y, y_err)`.
#[doc(alias = "gsl_fit_linear_est")]
pub fn linear_est(
    x: f64,
    c0: f64,
    c1: f64,
    cov00: f64,
    cov01: f64,
    cov11: f64,
) -> Result<(f64, f64), Error> {
    let mut y = 0.;
    let mut y_err = 0.;
    let ret =
        unsafe { sys::gsl_fit_linear_est(x, c0, c1, cov00, cov01, cov11, &mut y, &mut y_err) };
    Error::handle(ret, (y, y_err))
}

/// This function computes the best-fit linear regression coefficient c1 of the model Y = c_1 X for
/// the datasets (x, y), two vectors of length n with strides xstride and ystride.
/// The errors on y are assumed unknown so the variance of the parameter c1 is estimated from the
/// scatter of the points around the best-fit line and returned via the parameter cov11.
/// The sum of squares of the residuals from the best-fit line is returned in sumsq.
///
/// Returns `(c1, cov11, sumsq)`.
#[doc(alias = "gsl_fit_mul")]
pub fn mul<T: Vector<f64> + ?Sized>(x: &T, y: &T) -> Result<(f64, f64, f64), Error> {
    check_equal_len(x, y)?;
    let mut c1 = 0.;
    let mut cov11 = 0.;
    let mut sumsq = 0.;
    let ret = unsafe {
        sys::gsl_fit_mul(
            T::as_slice(x).as_ptr(),
            T::stride(x),
            T::as_slice(y).as_ptr(),
            T::stride(y),
            T::len(x),
            &mut c1,
            &mut cov11,
            &mut sumsq,
        )
    };
    Error::handle(ret, (c1, cov11, sumsq))
}

/// Returns `(c1, cov11, sumsq)`.
#[doc(alias = "gsl_fit_wmul")]
pub fn wmul<T: Vector<f64> + ?Sized>(x: &T, w: &T, y: &T) -> Result<(f64, f64, f64), Error> {
    check_equal_len(x, y)?;
    check_equal_len(x, w)?;
    let mut c1 = 0.;
    let mut cov11 = 0.;
    let mut sumsq = 0.;
    let ret = unsafe {
        sys::gsl_fit_wmul(
            T::as_slice(x).as_ptr(),
            T::stride(x),
            T::as_slice(w).as_ptr(),
            T::stride(w),
            T::as_slice(y).as_ptr(),
            T::stride(y),
            T::len(x),
            &mut c1,
            &mut cov11,
            &mut sumsq,
        )
    };
    Error::handle(ret, (c1, cov11, sumsq))
}

/// This function uses the best-fit linear regression coefficient c1 and its covariance cov11 to
/// compute the fitted function y and its standard deviation y_err for the model Y = c_1 X at the
/// point x.
///
/// Returns `(y, y_err)`.
#[doc(alias = "gsl_fit_mul_est")]
pub fn mul_est(x: f64, c1: f64, cov11: f64) -> Result<(f64, f64), Error> {
    let mut y = 0.;
    let mut y_err = 0.;
    let ret = unsafe { sys::gsl_fit_mul_est(x, c1, cov11, &mut y, &mut y_err) };
    Error::handle(ret, (y, y_err))
}
