//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# One dimensional Minimization

This module contains routines for finding minima of arbitrary
one-dimensional functions.  The library provides low level components
for a variety of iterative minimizers and convergence tests.  These
can be combined by the user to achieve the desired solution, with full
access to the intermediate steps of the algorithms. Each class of
methods uses the same framework, so that you can switch between
minimizers at runtime without needing to recompile your program. Each
instance of a minimizer keeps track of its own state, allowing the
minimizers to be used in multi-threaded programs.

## Overview

The minimization algorithms begin with a bounded region known to
contain a minimum.  The region is described by a lower bound $a$ and
an upper bound $b$, with an estimate of the location of the minimum
$x$ as shown in the next figure.

![Function with lower and upper bounds with an estimate of the minimum](https://www.gnu.org/software/gsl/doc/html/_images/min-interval.png)

The value of the function at $x$ must be less than the value of the
function at the ends of the interval,

$$f(a) > f(x) \text{ and } f(x) < f(b).$$

This condition guarantees that a minimum is contained somewhere within
the interval.  On each iteration a new point $x′$ is selected using
one of the available algorithms.  If the new point is a better
estimate of the minimum, i.e. where $f(x′) < f(x)$, then the current
estimate of the minimum $x$ is updated.  The new point also allows the
size of the bounded interval to be reduced, by choosing the most
compact set of points which satisfies the constraint $f(a) > f(x)
< f(b)$.  The interval is reduced until it encloses the true minimum
to a desired tolerance.  This provides a best estimate of the location
of the minimum and a rigorous error estimate.

Several bracketing algorithms are available within a single
framework. The user provides a high-level driver for the algorithm,
and the library provides the individual functions necessary for each
of the steps. There are three main phases of the iteration. The steps
are,

- initialize minimizer state, `s`, for algorithm `t`;
- update `s` using the iteration `t` with [`Minimizer::iterate`];
- test `s` for convergence (say using ), and repeat iteration if necessary.

The state for the minimizers is held in [`Minimizer`].  The updating
procedure uses only function evaluations (not derivatives).

## Caveats

Note that minimization functions can only search for one minimum at a
time.  When there are several minima in the search area, the first
minimum to be found will be returned; however it is difficult to
predict which of the minima this will be.  In most cases, *no error
will be reported if you try to find a minimum in an area where there
is more than one.*

With all minimization algorithms it can be difficult to determine the
location of the minimum to full numerical precision.  The behavior of
the function in the region of the minimum $x^*$ can be approximated by
a Taylor expansion,

$$y = f(x^\*) + (1/2) f′′(x^\*) (x - x^\*)²$$

and the second term of this expansion can be lost when added to the
first term at finite precision.  This magnifies the error in locating
$x^*$, making it proportional to $√ε$ (where $ε$ is the relative
accuracy of the floating point numbers).  For functions with higher
order minima, such as $x^4$, the magnification of the error is
correspondingly worse.  The best that can be achieved is to converge
to the limit of numerical accuracy in the function values, rather than
the location of the minimum itself.

## Stopping Parameters

A minimization procedure should stop when one of the following
conditions is true:

 * A minimum has been found to within the user-specified precision.
 * A user-specified maximum number of iterations has been reached.
 * An error has occurred.

The handling of these conditions is under user control. The function
[`test_interval`] allows the user to test the precision of the current
result.

## Minimization Algorithms

The minimization algorithms described in [`Type`] require an initial
interval which is guaranteed to contain a minimum—if $a$ and $b$ are
the endpoints of the interval and $x$ is an estimate of the minimum
then $f(a) > f(x) < f(b)$.  This ensures that the function has at
least one minimum somewhere in the interval.  If a valid initial
interval is used then these algorithm cannot fail, provided the
function is well-behaved.
*/

use crate::{Error, ffi::FFI, utilities::box_callback};
use std::ops::ControlFlow;

/// Type of the minimizer (without derivative).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Type {
    /// golden section algorithm.  See [`Minimizer::golden_section`].
    GoldenSection,
    /// Brent minimization algorithm which combines a parabolic
    /// interpolation with the golden section algorithm.  See
    /// [`Minimizer::brent`].
    Brent,
    /// variant of Brent’s algorithm which uses the safeguarded
    /// step-length algorithm of Gill and Murray.  See
    /// [`Minimizer::quad_golden`].
    QuadGolden,
}

impl Type {
    #[inline]
    fn to_c(self) -> *const sys::gsl_min_fminimizer_type {
        unsafe {
            match self {
                Self::GoldenSection => sys::gsl_min_fminimizer_goldensection,
                Self::Brent => sys::gsl_min_fminimizer_brent,
                Self::QuadGolden => sys::gsl_min_fminimizer_quad_golden,
            }
        }
    }
}

ffi_wrapper!(
    /// One dimensional minimization without derivatives.
    Minimizer<'a>,
    *mut sys::gsl_min_fminimizer,
    gsl_min_fminimizer_free
    ;f_struct: Option<Box<sys::gsl_function_struct>> => None;
    ;f: Option<Box<dyn Fn(f64) -> f64 + 'a>> => None;
);

impl<'a> Minimizer<'a> {
    /// This function returns a pointer to a newly allocated instance
    /// of a minimizer of type `t`.  For example, the following code
    /// creates an instance of a golden section minimizer,
    ///
    /// ```
    /// use rgsl::minimizer::{Minimizer, Type};
    /// let s = Minimizer::new(Type::GoldenSection);
    /// ```
    ///
    /// # Panic
    ///
    /// Panic if there is insufficient memory to create the minimizer.
    #[doc(alias = "gsl_min_fminimizer_alloc")]
    pub fn new(t: Type) -> Self {
        let ptr = unsafe { sys::gsl_min_fminimizer_alloc(t.to_c()) };

        if ptr.is_null() {
            panic!("rgsl::minimizer::Minimizer::new: out of memory")
        }
        Self::wrap(ptr)
    }

    /// The golden section algorithm is the simplest method of
    /// bracketing the minimum of a function.  It is the slowest
    /// algorithm provided by the library, with linear convergence.
    ///
    /// On each iteration, the algorithm first compares the
    /// sub-intervals from the endpoints to the current minimum.  The
    /// larger sub-interval is divided in a golden section (using the
    /// famous ratio $(3 - √5)/2 ≈ 0.3819660$ and the value of the
    /// function at this new point is calculated.  The new value is
    /// used with the constraint $f(a′) > f(x′) < f(b′)$ to a select
    /// new interval containing the minimum, by discarding the least
    /// useful point.  This procedure can be continued indefinitely
    /// until the interval is sufficiently small.  Choosing the golden
    /// section as the bisection ratio can be shown to provide the
    /// fastest convergence for this type of algorithm.
    #[doc(alias = "gsl_min_fminimizer_goldensection")]
    pub fn golden_section() -> Self {
        Self::new(Type::GoldenSection)
    }

    /// The Brent minimization algorithm combines a parabolic
    /// interpolation with the golden section algorithm.  This
    /// produces a fast algorithm which is still robust.
    ///
    /// The outline of the algorithm can be summarized as follows: on
    /// each iteration Brent’s method approximates the function using
    /// an interpolating parabola through three existing points.  The
    /// minimum of the parabola is taken as a guess for the minimum.
    /// If it lies within the bounds of the current interval then the
    /// interpolating point is accepted, and used to generate a
    /// smaller interval.  If the interpolating point is not accepted
    /// then the algorithm falls back to an ordinary golden section
    /// step.  The full details of Brent’s method include some
    /// additional checks to improve convergence.
    #[doc(alias = "gsl_min_fminimizer_brent")]
    pub fn brent() -> Self {
        Self::new(Type::Brent)
    }

    /// This is a variant of Brent’s algorithm which uses the
    /// safeguarded step-length algorithm of Gill and Murray.
    #[doc(alias = "gsl_min_fminimizer_quad_golden")]
    pub fn quad_golden() -> Self {
        Self::new(Type::QuadGolden)
    }

    /// This function sets, or resets, an existing minimizer s to use
    /// the function f and the initial search interval [x_lower,
    /// x_upper], with a guess for the location of the minimum
    /// x_minimum.
    ///
    /// If the interval given does not contain a minimum, then the
    /// function returns an error code of `Error::Invalid`.
    #[doc(alias = "gsl_min_fminimizer_set")]
    pub fn set<F: Fn(f64) -> f64 + 'a>(
        &mut self,
        f: F,
        x_minimum: f64,
        x_lower: f64,
        x_upper: f64,
    ) -> Result<(), Error> {
        let mut f = Box::new(f);
        let mut f_struct = unsafe { box_callback(&mut f) };
        self.f = Some(f);

        let ret = unsafe {
            sys::gsl_min_fminimizer_set(
                self.unwrap_unique(),
                &mut *f_struct,
                x_minimum,
                x_lower,
                x_upper,
            )
        };
        self.f_struct = Some(f_struct);
        Error::handle(ret, ())
    }

    /// This function is equivalent to gsl_min_fminimizer_set but uses the values f_minimum, f_lower
    /// and f_upper instead of computing f(x_minimum), f(x_lower) and f(x_upper).
    #[doc(alias = "gsl_min_fminimizer_set_with_values")]
    pub fn set_with_values<F: Fn(f64) -> f64 + 'a>(
        &mut self,
        f: F,
        x_minimum: f64,
        f_minimum: f64,
        x_lower: f64,
        f_lower: f64,
        x_upper: f64,
        f_upper: f64,
    ) -> Result<(), Error> {
        let mut f = Box::new(f);
        let mut f_struct = unsafe { box_callback(&mut f) };
        self.f = Some(f);

        let ret = unsafe {
            sys::gsl_min_fminimizer_set_with_values(
                self.unwrap_unique(),
                &mut *f_struct,
                x_minimum,
                f_minimum,
                x_lower,
                f_lower,
                x_upper,
                f_upper,
            )
        };
        self.f_struct = Some(f_struct);
        Error::handle(ret, ())
    }

    #[doc(alias = "gsl_min_fminimizer_name")]
    pub fn name(&self) -> Type {
        let n = unsafe { sys::gsl_min_fminimizer_name(self.unwrap_shared()) };
        map_name!(
            rgsl::minimizer::Minimizer,
            [
                (c"goldensection", Type::GoldenSection),
                (c"brent", Type::Brent),
                (c"quad-golden", Type::QuadGolden)
            ],
            n,
            Type
        )
    }

    /// Return the current estimate of the position of the minimum.
    #[doc(alias = "gsl_min_fminimizer_x_minimum")]
    pub fn x_minimum(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_x_minimum(self.unwrap_shared()) }
    }

    /// Return the current lower bound of the interval where the minimum lies.
    #[doc(alias = "gsl_min_fminimizer_x_lower")]
    pub fn x_lower(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_x_lower(self.unwrap_shared()) }
    }

    /// Return the current upper bound of the interval where the minimum lies.
    #[doc(alias = "gsl_min_fminimizer_x_upper")]
    pub fn x_upper(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_x_upper(self.unwrap_shared()) }
    }

    /// Return the value of the function at the current estimate of
    /// the minimum of the interval.
    #[doc(alias = "gsl_min_fminimizer_f_minimum")]
    pub fn f_minimum(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_f_minimum(self.unwrap_shared()) }
    }

    /// Return the value of the function at the lower bound of the interval.
    #[doc(alias = "gsl_min_fminimizer_f_lower")]
    pub fn f_lower(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_f_lower(self.unwrap_shared()) }
    }

    /// Return the value of the function at the upper bound of the interval.
    #[doc(alias = "gsl_min_fminimizer_f_upper")]
    pub fn f_upper(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_f_upper(self.unwrap_shared()) }
    }

    #[deprecated(since = "8.0.0", note = "Use x_minimum")]
    #[doc(alias = "gsl_min_fminimizer_minimum")]
    pub fn minimum(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_minimum(self.unwrap_shared()) }
    }

    /// This function performs a single iteration of the minimizer
    /// s. If the iteration encounters an unexpected problem then an
    /// error code will be returned,
    ///
    /// - [`Error::BadFunction`] the iteration encountered a singular
    ///   point where the function evaluated to Inf or NaN.
    /// - [`Error::Failure`] the algorithm could not improve the
    ///   current best approximation or bounding interval.
    ///
    /// The minimizer maintains a current best estimate of the
    /// position of the minimum at all times, and the current interval
    /// bounding the minimum. This information can be accessed with
    /// the following auxiliary functions,
    #[doc(alias = "gsl_min_fminimizer_iterate")]
    pub fn iterate(&mut self) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_min_fminimizer_iterate(self.unwrap_unique()) };
        Error::handle(ret, ())
    }
}

/// This function tests for the convergence of the interval
/// \[`x_lower`, `x_upper`\] with absolute error epsabs and relative
/// error epsrel. The test returns `ControlFlow::Break(())` if the
/// following condition is achieved,
///
/// $$|a - b| < \epsabs + \epsrel \min(|a|,|b|)$$
///
/// when the interval $x = \[a,b\]$ does not include the origin.  If
/// the interval includes the origin then $\min(|a|,|b|)$ is replaced
/// by zero (which is the minimum value of |x| over the interval).
/// This ensures that the relative error is accurately estimated for
/// minima close to the origin.
///
/// This condition on the interval also implies that any estimate of
/// the minimum $x_m$ in the interval satisfies the same condition with
/// respect to the true minimum $x_m^*$,
///
/// $$|x_m - x_m^*| < epsabs + epsrel x_m^*$$
///
/// assuming that the true minimum $x_m^*$ is contained within the interval.
#[doc(alias = "gsl_min_test_interval")]
pub fn test_interval(x_lower: f64, x_upper: f64, epsabs: f64, epsrel: f64) -> ControlFlow<()> {
    Error::control_flow(unsafe { sys::gsl_min_test_interval(x_lower, x_upper, epsabs, epsrel) })
}

#[cfg(any(test, doctest))]
mod test {
    /// This doc block will be used to ensure that the closure can't
    /// be set everywhere!
    ///
    /// ```compile_fail
    /// use crate::rgsl::*;
    /// use crate::rgsl::minimizer::test_interval;
    ///
    /// fn set(min: &mut Minimizer) {
    ///     let y = "lalal".to_owned();
    ///     min.set(|x| {
    ///         println!("==> {:?}", y);
    ///         x * x - 5.
    ///     }, 1.0, -5.0, 5.0);
    /// }
    ///
    /// let mut min = Minimizer::brent();
    /// set(&mut min);
    /// let status = min.iterate();
    /// ```
    ///
    /// Same but a working version:
    ///
    /// ```
    /// use crate::rgsl::*;
    /// use crate::rgsl::minimizer::test_interval;
    ///
    /// fn set(min: &mut Minimizer) {
    ///     min.set(|x| x * x - 5., 1.0, -5.0, 5.0);
    /// }
    ///
    /// let mut min = Minimizer::brent();
    /// set(&mut min);
    /// let status = min.iterate();
    /// ```
    use super::*;
    use crate::minimizer::test_interval;

    #[test]
    fn test_minizer_pointers() -> Result<(), Error> {
        // This function will move the `Minimizer`, testing that its
        // inner pointers stay valid.
        fn create() -> Result<Minimizer<'static>, Error> {
            let mut m = Minimizer::golden_section();
            m.set(|x| x.powi(2), 0., -1., 1.)?;
            Ok(m)
        }
        let mut m = create()?;
        m.iterate()?;
        Ok(())
    }

    #[test]
    fn test_names() {
        let m = Minimizer::golden_section();
        assert_eq!(m.name(), Type::GoldenSection);
        let min = Minimizer::brent();
        assert_eq!(min.name(), Type::Brent);
        let min = Minimizer::quad_golden();
        assert_eq!(min.name(), Type::QuadGolden);
    }

    fn quadratic_test_fn(x: f64) -> f64 {
        x.powf(2.0) - 5.0
    }

    #[test]
    fn test_min() -> Result<(), Error> {
        let mut min = Minimizer::brent();
        min.set(quadratic_test_fn, 1.0, -5.0, 5.0).unwrap();

        let max_iter = 5_usize;
        let eps_abs = 0.0001;
        let eps_rel = 0.0000001;

        for iter in 0..max_iter {
            // iterate for next value
            min.iterate()?;

            // test for convergence
            let r = min.x_minimum();
            let x_lo = min.x_lower();
            let x_hi = min.x_upper();

            // print current iteration
            println!("{} [{}, {}] {} {}", iter, x_lo, x_hi, r, x_hi - x_lo);

            let status = test_interval(x_lo, x_hi, eps_abs, eps_rel);
            if matches!(status, ControlFlow::Break(())) {
                println!("Converged");
                break;
            }
        }
        Ok(())
    }
}
