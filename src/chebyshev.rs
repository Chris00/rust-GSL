//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Chebyshev Approximations

This module describes routines for computing Chebyshev approximations
to univariate functions.  A Chebyshev approximation is a truncation of
the series $f(x) = ∑ cₙ Tₙn(x)$, where the Chebyshev polynomials
$Tₙ(x) = \cos(n \arccos x)$ provide an orthogonal basis of polynomials
on the interval $[-1,1]$ with the weight function $1 / √{1-x^2}$.  The
first few Chebyshev polynomials are, $T₀(x) = 1$, $T₁(x) = x$, $T₂(x)
= 2 x^2 - 1$.

For further information see Abramowitz & Stegun, Chapter 22.

## References and Further Reading

The following paper describes the use of Chebyshev series,

R. Broucke, “Ten Subroutines for the Manipulation of Chebyshev Series
\[C1\] (Algorithm 446)”.  Communications of the ACM 16(4), 254–256
(1973).
*/

use crate::Error;
use crate::{ffi::FFI, utilities::wrap_callback};

ffi_wrapper!(
    /// Storage for a Chebyshev series.
    ///
    /// # Example
    ///
    /// ```
    /// use rgsl::chebyshev::ChebSeries;
    /// let f = |x: f64| if x < 0.5 { 0.25 } else { 0.75 };
    /// let mut cs = ChebSeries::new(40);
    /// cs.init(f, 0., 1.)?;
    /// let n = 10000;
    /// for i in 0..n {
    ///     let x = i as f64 / n as f64;
    ///     let r10 = cs.eval_n(10, x);
    ///     let r40 = cs.eval(x);
    ///     println!("{x} {} {r10} {r40}", f(x));
    /// }
    /// # Ok::<(), rgsl::Error>(())
    /// ```
    ChebSeries,
    *mut sys::gsl_cheb_series, gsl_cheb_free);

impl ChebSeries {
    /// Return a non-initialized space for a Chebyshev series of order
    /// `n`.
    ///
    /// # Panic
    /// Panic if there is not enough memory.
    #[doc(alias = "gsl_cheb_alloc")]
    pub fn new(n: usize) -> Self {
        let tmp = unsafe { sys::gsl_cheb_alloc(n) };

        if tmp.is_null() {
            panic!("rgsl::chebyshev::ChebSeries::new: out of memory")
        }
        Self::wrap(tmp)
    }

    /// Compute the Chebyshev approximation cs for the function `f`
    /// over the range (`a`, `b`) to the order specified in
    /// [`Self::new`].  The computation of the Chebyshev approximation
    /// is an $O(n²)$ process, and requires $n$ function evaluations.
    ///
    /// The approximation is made over the range $\[a,b\]$ using `order + 1`
    /// terms, including the coefficient $c₀$.  The series is computed using
    /// the following convention,
    ///
    /// $$f(x) = (c₀ / 2) + ∑_{n=1}^∞ cₙ Tₙ(x)$$
    ///
    /// which is needed when accessing the coefficients directly.
    ///
    /// # Panic
    /// Panic if `a >= b`.
    #[doc(alias = "gsl_cheb_init")]
    pub fn init<F: FnMut(f64) -> f64>(&mut self, mut f: F, a: f64, b: f64) -> Result<(), Error> {
        // The `function` is not stored in `Self` so it is fine to
        // keep it on the stack.
        let function = unsafe { wrap_callback(&mut f) };

        let ret = unsafe { sys::gsl_cheb_init(self.unwrap_unique(), &function, a, b) };
        Error::handle(ret, ())
    }

    /// Return the order of Chebyshev series.
    #[doc(alias = "gsl_cheb_order")]
    pub fn order(&self) -> usize {
        unsafe { sys::gsl_cheb_order(self.unwrap_shared()) }
    }

    /// Return the size of the Chebyshev coefficient array for the
    /// Chebyshev series.
    #[doc(alias = "gsl_cheb_size")]
    pub fn size(&self) -> usize {
        unsafe { sys::gsl_cheb_size(self.unwrap_shared()) }
    }

    /// Evaluate the Chebyshev series at a given point `x`.
    #[doc(alias = "gsl_cheb_eval")]
    pub fn eval(&self, x: f64) -> f64 {
        unsafe { sys::gsl_cheb_eval(self.unwrap_shared(), x) }
    }

    /// Compute the Chebyshev series at a given point `x`, returning
    /// `(r, abserr)` where `r` is the series result and `abserr` is
    /// its absolute error.  The error estimate is made from the first
    /// neglected term in the series.
    #[doc(alias = "gsl_cheb_eval_err")]
    pub fn eval_err(&self, x: f64) -> Result<(f64, f64), Error> {
        let mut result = 0.;
        let mut abs_err = 0.;

        let ret =
            unsafe { sys::gsl_cheb_eval_err(self.unwrap_shared(), x, &mut result, &mut abs_err) };
        Error::handle(ret, (result, abs_err))
    }

    /// Evaluate the Chebyshev series at a given point `x`, to (at
    /// most) the given order `order`.
    #[doc(alias = "gsl_cheb_eval_n")]
    pub fn eval_n(&self, order: usize, x: f64) -> f64 {
        unsafe { sys::gsl_cheb_eval_n(self.unwrap_shared(), order, x) }
    }

    /// Evaluate the Chebyshev series at a given point `x`, returning
    /// `(r, abserr)` where `r` is the series result and `abserr` is
    /// its absolute error, to (at most) the given order `order`.  The
    /// error estimate is made from the first neglected term in the
    /// series.
    #[doc(alias = "gsl_cheb_eval_n_err")]
    pub fn eval_n_err(&self, order: usize, x: f64) -> Result<(f64, f64), Error> {
        let mut result = 0.;
        let mut abs_err = 0.;

        let ret = unsafe {
            sys::gsl_cheb_eval_n_err(self.unwrap_shared(), order, x, &mut result, &mut abs_err)
        };
        Error::handle(ret, (result, abs_err))
    }

    /// Compute the derivative of the Chebyshev series, storing the
    /// derivative coefficients in the `deriv`.  The two series `self`
    /// and `deriv` must have been allocated with the same order.
    #[doc(alias = "gsl_cheb_calc_deriv")]
    pub fn calc_deriv(&self, deriv: &mut ChebSeries) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_cheb_calc_deriv(deriv.unwrap_unique(), self.unwrap_shared()) };
        Error::handle(ret, ())
    }

    /// Compute the integral of the Chebyshev series, storing the
    /// integral coefficients in `integ`.  The two series `self` and
    /// `integ` must have been allocated with the same order.  The
    /// lower limit of the integration is taken to be the left hand
    /// end of the range `a` (see [`Self::init`]).
    #[doc(alias = "gsl_cheb_calc_integ")]
    pub fn calc_integ(&self, integ: &mut ChebSeries) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_cheb_calc_integ(integ.unwrap_unique(), self.unwrap_shared()) };
        Error::handle(ret, ())
    }
}
