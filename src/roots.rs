//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# One dimensional Root-Finding

This module contains routines for finding roots of arbitrary
one-dimensional functions.  The library provides low level components
for a variety of iterative solvers and convergence tests.  These can
be combined by the user to achieve the desired solution, with full
access to the intermediate steps of the iteration.  Each class of
methods uses the same framework, so that you can switch between
solvers at runtime without needing to recompile your program.  Each
instance of a solver keeps track of its own state, allowing the
solvers to be used in multi-threaded programs.

## Overview

One-dimensional root finding algorithms can be divided into two
classes, *root bracketing* and *root polishing*.  Algorithms which
proceed by bracketing a root are guaranteed to converge.  Bracketing
algorithms begin with a bounded region known to contain a root.  The
size of this bounded region is reduced, iteratively, until it encloses
the root to a desired tolerance.  This provides a rigorous error
estimate for the location of the root.

The technique of root polishing attempts to improve an initial guess
to the root.  These algorithms converge only if started “close enough”
to a root, and sacrifice a rigorous error bound for speed.  By
approximating the behavior of a function in the vicinity of a root
they attempt to find a higher order improvement of an initial guess.
When the behavior of the function is compatible with the algorithm and
a good initial guess is available a polishing algorithm can provide
rapid convergence.

In GSL both types of algorithm are available in similar frameworks.
The user provides a high-level driver for the algorithms, and the
library provides the individual functions necessary for each of the
steps.  There are three main phases of the iteration. The steps are,
- initialize solver state, `s`, for algorithm `t`;
- update `s` using the iteration `t`;
- test `s` for convergence, and repeat iteration if necessary.

The state for bracketing solvers is held in a [`RootFSolver`].  The
updating procedure uses only function evaluations (not derivatives).
The state for root polishing solvers is held in a [`RootFdfSolver`].
The updates require both the function and its derivative (hence the
name “fdf”) to be supplied by the user.

# Search Stopping Parameters

A root finding procedure should stop when one of the following
conditions is true:

- A root has been found to within the user-specified precision.
- A user-specified maximum number of iterations has been reached.
- An error has occurred.

The handling of these conditions is under user control.  The functions
[`test_interval`], [`test_delta`], and [`test_residual`] allow the
user to test the precision of the current result in several standard
ways.
 */

use crate::{Error, ffi::FFI, utilities::box_callback};
use std::{
    ffi::{c_double, c_void},
    ops::ControlFlow,
};

/// Type of root bracketing algorithm.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Type {
    /// Bisection algorithm.  See [`RootFSolver::bisection`].
    Bisection,
    /// False position algorithm.  See [`RootFSolver::falsepos`].
    FalsePos,
    /// Brent-Dekker method.  See [`RootFSolver::brent`].
    Brent,
}

impl Type {
    #[inline]
    fn to_c(self) -> *const sys::gsl_root_fsolver_type {
        unsafe {
            match self {
                Self::Bisection => sys::gsl_root_fsolver_bisection,
                Self::FalsePos => sys::gsl_root_fsolver_falsepos,
                Self::Brent => sys::gsl_root_fsolver_brent,
            }
        }
    }
}

ffi_wrapper!(
    /// Root bracketing algorithms solver.
    ///
    /// Root bracketing algorithms require an initial interval which
    /// is guaranteed to contain a root—if $a$ and $b$ are the
    /// endpoints of the interval then $f(a)$ must differ in sign from
    /// $f(b)$.  This ensures that the function crosses zero at least
    /// once in the interval.  If a valid initial interval is used
    /// then these algorithm cannot fail, provided the function is
    /// well-behaved.  Note that a bracketing algorithm cannot find
    /// roots of even degree, since these do not cross the x-axis.
    RootFSolver<'a>,
    *mut sys::gsl_root_fsolver,
    gsl_root_fsolver_free
    ;f_struct: Option<Box<sys::gsl_function_struct>> => None;
    ;f: Option<Box<dyn FnMut(f64) -> f64 + 'a>> => None;
);

impl<'a> RootFSolver<'a> {
    /// Allocate an instance of a solver of type `t`.
    ///
    /// # Panic
    /// Panic If there is insufficient memory to create the solver.
    #[doc(alias = "gsl_root_fsolver_alloc")]
    pub fn new(t: Type) -> Self {
        let tmp = unsafe { sys::gsl_root_fsolver_alloc(t.to_c()) };

        if tmp.is_null() {
            panic!("rgsl::roots::RootFSolver::new: out of memory");
        }
        Self::wrap(tmp)
    }

    /// The bisection algorithm is the simplest method of bracketing
    /// the roots of a function.  It is the slowest algorithm provided
    /// by the library, with linear convergence.  On each iteration,
    /// the interval is bisected and the value of the function at the
    /// midpoint is calculated. The sign of this value is used to
    /// determine which half of the interval does not contain a
    /// root. That half is discarded to give a new, smaller interval
    /// containing the root. This procedure can be continued
    /// indefinitely until the interval is sufficiently small.
    ///
    /// At any time the current estimate of the root is taken as the
    /// midpoint of the interval.
    #[doc(alias = "gsl_root_fsolver_bisection")]
    pub fn bisection() -> Self {
        Self::new(Type::Bisection)
    }

    /// The false position algorithm is a method of finding roots
    /// based on linear interpolation.  Its convergence is linear, but
    /// it is usually faster than bisection.
    ///
    /// On each iteration a line is drawn between the endpoints $(a,
    /// f(a))$ and $(b, f(b))$ and the point where this line crosses
    /// the x-axis taken as a “midpoint”.  The value of the function
    /// at this point is calculated and its sign is used to determine
    /// which side of the interval does not contain a root.  That side
    /// is discarded to give a new, smaller interval containing the
    /// root.  This procedure can be continued indefinitely until the
    /// interval is sufficiently small.
    ///
    /// The best estimate of the root is taken from the linear
    /// interpolation of the interval on the current iteration.
    #[doc(alias = "gsl_root_fsolver_falsepos")]
    pub fn falsepos() -> Self {
        Self::new(Type::FalsePos)
    }

    /// The Brent-Dekker method (referred to here as Brent’s method)
    /// combines an interpolation strategy with the bisection
    /// algorithm.  This produces a fast algorithm which is still
    /// robust.
    ///
    /// On each iteration Brent’s method approximates the function
    /// using an interpolating curve.  On the first iteration this is
    /// a linear interpolation of the two endpoints.  For subsequent
    /// iterations the algorithm uses an inverse quadratic fit to the
    /// last three points, for higher accuracy.  The intercept of the
    /// interpolating curve with the x-axis is taken as a guess for
    /// the root.  If it lies within the bounds of the current
    /// interval then the interpolating point is accepted, and used to
    /// generate a smaller interval.  If the interpolating point is
    /// not accepted then the algorithm falls back to an ordinary
    /// bisection step.
    ///
    /// The best estimate of the root is taken from the most recent
    /// interpolation or bisection.
    #[doc(alias = "gsl_root_fsolver_brent")]
    pub fn brent() -> Self {
        Self::new(Type::Brent)
    }

    /// Initialize, or reinitialize, the solver to use the function `f`
    /// and the initial search interval `\[x_lower, x_upper\]`.
    #[doc(alias = "gsl_root_fsolver_set")]
    pub fn set<F: FnMut(f64) -> f64 + 'a>(
        &mut self,
        f: F,
        x_lower: f64,
        x_upper: f64,
    ) -> Result<(), Error> {
        let mut f = Box::new(f);
        let mut f_struct = unsafe { box_callback(&mut f) };
        self.f = Some(f);

        let ret = unsafe {
            sys::gsl_root_fsolver_set(self.unwrap_unique(), &mut *f_struct, x_lower, x_upper)
        };
        self.f_struct = Some(f_struct);
        Error::handle(ret, ())
    }

    /// Performs one iteration of the solver to update its state
    /// (according to the solver's type).  If the iteration encounters
    /// an unexpected problem then an error code will be returned.
    ///
    /// The solver maintains a current best estimate of the root at
    /// all times.  The bracketing solvers also keep track of the
    /// current best interval bounding the root.
    #[doc(alias = "gsl_root_fsolver_iterate")]
    pub fn iterate(&mut self) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_root_fsolver_iterate(self.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// Returns the solver type name.
    #[doc(alias = "gsl_root_fsolver_name")]
    pub fn name(&self) -> Type {
        let n = unsafe { sys::gsl_root_fsolver_name(self.unwrap_shared()) };
        map_name!(
            rgsl::roots::RootFSolver,
            [
                (c"bisection", Type::Bisection),
                (c"falsepos", Type::FalsePos),
                (c"brent", Type::Brent),
            ],
            n,
            Type
        )
    }

    /// Return the current estimate of the root.
    #[doc(alias = "gsl_root_fsolver_root")]
    pub fn root(&self) -> f64 {
        unsafe { sys::gsl_root_fsolver_root(self.unwrap_shared()) }
    }

    /// Return the current lower bound of the bracketing interval.
    #[doc(alias = "gsl_root_fsolver_x_lower")]
    pub fn x_lower(&self) -> f64 {
        unsafe { sys::gsl_root_fsolver_x_lower(self.unwrap_shared()) }
    }

    /// Return the current upper bound of the bracketing interval.
    #[doc(alias = "gsl_root_fsolver_x_upper")]
    pub fn x_upper(&self) -> f64 {
        unsafe { sys::gsl_root_fsolver_x_upper(self.unwrap_shared()) }
    }
}

/// Type of root polishing algorithm.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TypeFdf {
    /// Newton algorithm.  See [`RootFdfSolver::newton`].
    Newton,
    /// Secant method.  See [`RootFdfSolver::secant`].
    Secant,
    /// Steffenen method.  See [`RootFdfSolver::steffensen`].
    Steffensen,
}

impl TypeFdf {
    fn to_c(self) -> *const sys::gsl_root_fdfsolver_type {
        unsafe {
            match self {
                Self::Newton => sys::gsl_root_fdfsolver_newton,
                Self::Secant => sys::gsl_root_fdfsolver_secant,
                Self::Steffensen => sys::gsl_root_fdfsolver_steffenson,
            }
        }
    }
}

/// Specify which data-types may be used to specify a function and its
/// derivatives.
pub trait Fdf {
    /// Return the value of the function $f$ at `x`.
    fn f(&mut self, x: f64) -> f64;
    /// Return the value of the function $f′$ at `x`.
    fn df(&mut self, x: f64) -> f64;
    /// Return the both value of the function $f$ and $f′$ at `x`.
    fn fdf(&mut self, x: f64) -> (f64, f64);
}

impl<F: FnMut(f64) -> f64, DF: FnMut(f64) -> f64> Fdf for (F, DF) {
    fn f(&mut self, x: f64) -> f64 {
        self.0(x)
    }
    fn df(&mut self, x: f64) -> f64 {
        self.1(x)
    }
    fn fdf(&mut self, x: f64) -> (f64, f64) {
        (self.0(x), self.1(x))
    }
}

impl<FDF: FnMut(f64) -> (f64, f64)> Fdf for FDF {
    fn f(&mut self, x: f64) -> f64 {
        self(x).0
    }
    fn df(&mut self, x: f64) -> f64 {
        self(x).1
    }
    fn fdf(&mut self, x: f64) -> (f64, f64) {
        self(x)
    }
}

impl<F, DF, FDF> Fdf for (F, DF, FDF)
where
    F: FnMut(f64) -> f64,
    DF: FnMut(f64) -> f64,
    FDF: FnMut(f64) -> (f64, f64),
{
    fn f(&mut self, x: f64) -> f64 {
        self.0(x)
    }
    fn df(&mut self, x: f64) -> f64 {
        self.1(x)
    }
    fn fdf(&mut self, x: f64) -> (f64, f64) {
        self.2(x)
    }
}

ffi_wrapper!(
    /// Root finding methods which do require derivatives.
    ///
    /// The root polishing algorithms require an initial guess for the
    /// location of the root.  There is no absolute guarantee of
    /// convergence—the function must be suitable for this technique
    /// and the initial guess must be sufficiently close to the root
    /// for it
    ///
    /// When these conditions are satisfied then convergence is
    /// quadratic.  These algorithms make use of both the function and
    /// its derivative.
    RootFdfSolver<'a>,
    *mut sys::gsl_root_fdfsolver,
    gsl_root_fdfsolver_free
    ;fdf_struct: Option<Box<sys::gsl_function_fdf_struct>> => None;
    ;fdf: Option<Box<dyn Fdf + 'a>> => None;);

impl<'a> RootFdfSolver<'a> {
    /// Return a new instance of a derivative-based solver of type `t`.
    ///
    /// # Panic
    /// Panic if there is insufficient memory to create the solver.
    #[doc(alias = "gsl_root_fdfsolver_alloc")]
    pub fn new(t: TypeFdf) -> RootFdfSolver<'a> {
        let tmp = unsafe { sys::gsl_root_fdfsolver_alloc(t.to_c()) };

        if tmp.is_null() {
            panic!("rgsl::roots::RootFdfSolver::new: out of memory")
        }
        Self::wrap(tmp)
    }

    /// Newton’s Method is the standard root-polishing algorithm.  The
    /// algorithm begins with an initial guess for the location of the
    /// root. On each iteration, a line tangent to the function f is
    /// drawn at that position. The point where this line crosses the
    /// x-axis becomes the new guess.  The iteration is defined by the
    /// following sequence,
    ///
    /// $$x_{i+1} = x_i - \frac{f(x_i)}{f′(x_i)}$$
    ///
    /// Newton’s method converges quadratically for single roots, and
    /// linearly for multiple roots.
    #[doc(alias = "gsl_root_fdfsolver_newton")]
    pub fn newton() -> Self {
        Self::new(TypeFdf::Newton)
    }

    /// The secant method is a simplified version of Newton’s method
    /// which does not require the computation of the derivative on
    /// every step.
    ///
    /// On its first iteration the algorithm begins with Newton’s
    /// method, using the derivative to compute a first step,
    ///
    /// $$x_1 = x_0 - \frac{f(x_0)}{f′(x_0)}$$
    ///
    /// Subsequent iterations avoid the evaluation of the derivative
    /// by replacing it with a numerical estimate, the slope of the
    /// line through the previous two points,
    ///
    /// $$x_{i+1} = x_i - \frac{f(x_i)}{f′_{est}}
    /// \quad\text{where }
    /// f′_{est} = \frac{f(x_{i}) - f(x_{i-1})}{x_i - x_{i-1}}$$
    ///
    /// When the derivative does not change significantly in the
    /// vicinity of the root the secant method gives a useful saving.
    /// Asymptotically the secant method is faster than Newton’s
    /// method whenever the cost of evaluating the derivative is more
    /// than 0.44 times the cost of evaluating the function itself.
    /// As with all methods of computing a numerical derivative the
    /// estimate can suffer from cancellation errors if the separation
    /// of the points becomes too small.
    ///
    /// On single roots, the method has a convergence of order (1 +
    /// √5)/2 (approximately 1.62).  It converges linearly for
    /// multiple roots.
    #[doc(alias = "gsl_root_fdfsolver_secant")]
    pub fn secant() -> Self {
        Self::new(TypeFdf::Secant)
    }

    /// The Steffensen Method¹ provides the fastest convergence of all
    /// the routines.  It combines the basic Newton algorithm with an
    /// Aitken “delta-squared” acceleration.  If the Newton iterates
    /// are $x_i$ then the acceleration procedure generates a new
    /// sequence $R_i$,
    ///
    /// $$R_i = x_i - \frac{(x_{i+1} - x_i)²}{(x_{i+2} - 2 x_{i+1} + x_i)}$$
    ///
    /// which converges faster than the original sequence under
    /// reasonable conditions.  The new sequence requires three terms
    /// before it can produce its first value so the method returns
    /// accelerated values on the second and subsequent iterations.
    /// On the first iteration it returns the ordinary Newton
    /// estimate. The Newton iterate is also returned if the
    /// denominator of the acceleration term ever becomes zero.
    ///
    /// As with all acceleration procedures this method can become
    /// unstable if the function is not well-behaved.
    ///
    /// ¹ J. F. Steffensen (1873–1961). The spelling used in the name
    /// of the function is slightly incorrect, but has been preserved
    /// to avoid incompatibility.
    #[doc(alias = "gsl_root_fdfsolver_steffenson")]
    pub fn steffensen() -> Self {
        Self::new(TypeFdf::Steffensen)
    }

    /// This function initializes, or reinitializes, the solver to use
    /// the function and derivative fdf and the initial guess root.
    #[doc(alias = "gsl_root_fdfsolver_set")]
    pub fn set<FDF: Fdf + 'a>(&mut self, fdf: FDF, root: f64) -> Result<(), Error> {
        unsafe extern "C" fn inner_f<FDF: Fdf>(x: c_double, params: *mut c_void) -> f64 {
            let t = unsafe { &mut *params.cast::<FDF>() };
            FDF::f(t, x)
        }

        unsafe extern "C" fn inner_df<FDF: Fdf>(x: c_double, params: *mut c_void) -> f64 {
            let t = unsafe { &mut *params.cast::<FDF>() };
            FDF::df(t, x)
        }

        unsafe extern "C" fn inner_fdf<FDF: Fdf>(
            x: c_double,
            params: *mut c_void,
            y: *mut c_double,
            dy: *mut c_double,
        ) {
            let t = unsafe { &mut *params.cast::<FDF>() };
            let (fx, dfx) = FDF::fdf(t, x);
            unsafe {
                *y = fx;
                *dy = dfx;
            }
        }

        let mut fdf = Box::new(fdf);
        let mut fdf_struct = Box::new(sys::gsl_function_fdf {
            f: Some(inner_f::<FDF>),
            df: Some(inner_df::<FDF>),
            fdf: Some(inner_fdf::<FDF>),
            params: &mut *fdf as *mut FDF as *mut _,
        });
        self.fdf = Some(fdf);

        let ret =
            unsafe { sys::gsl_root_fdfsolver_set(self.unwrap_unique(), &mut *fdf_struct, root) };
        self.fdf_struct = Some(fdf_struct);
        Error::handle(ret, ())
    }

    /// Perform a single iteration of the solver (according to its
    /// type).  If the iteration encounters an unexpected problem then
    /// an error code will be returned.
    ///
    /// The solver maintains a current best estimate of the root at
    /// all times.  The bracketing solvers also keep track of the
    /// current best interval bounding the root.
    #[doc(alias = "gsl_root_fdfsolver_iterate")]
    pub fn iterate(&mut self) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_root_fdfsolver_iterate(self.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// Return the solver type name.
    #[doc(alias = "gsl_root_fdfsolver_name")]
    pub fn name(&self) -> TypeFdf {
        let n = unsafe { sys::gsl_root_fdfsolver_name(self.unwrap_shared()) };
        map_name!(
            rgsl::roots::RootFdfSolver,
            [
                (c"newton", TypeFdf::Newton),
                (c"secant", TypeFdf::Secant),
                (c"steffenson", TypeFdf::Steffensen),
            ],
            n,
            TypeFdf
        )
    }

    /// Return the current estimate of the root.
    #[doc(alias = "gsl_root_fdfsolver_root")]
    pub fn root(&self) -> f64 {
        unsafe { sys::gsl_root_fdfsolver_root(self.unwrap_shared()) }
    }
}

/// Test for the convergence of the interval \[`x_lower`, `x_upper`\]
/// with absolute error `epsabs` and relative error `epsrel`.
///
/// The test returns [`ControlFlow::Break`] if the following condition
/// is achieved,
///
/// $$|a - b| < \epsabs + \epsrel \\, \min(|a|,|b|)$$
///
/// when the interval $x = \[a,b\]$ does not include the origin.  If
/// the interval includes the origin then $\min(|a|,|b|)$ is replaced
/// by zero (which is the minimum value of $|x|$ over the interval).
/// This ensures that the relative error is accurately estimated for
/// roots close to the origin.
///
/// This condition on the interval also implies that any estimate of
/// the root $r$ in the interval satisfies the same condition with
/// respect to the true root $r^*$,
///
/// $$|r - r^\*| < \epsabs + \epsrel \\, r^\*$$
///
/// assuming that the true root r^* is contained within the interval.
#[doc(alias = "gsl_root_test_interval")]
pub fn test_interval(x_lower: f64, x_upper: f64, epsabs: f64, epsrel: f64) -> ControlFlow<()> {
    Error::control_flow(unsafe { sys::gsl_root_test_interval(x_lower, x_upper, epsabs, epsrel) })
}

/// Test for the convergence of the sequence `x0`, `x1` with absolute
/// error `epsabs` and relative error `epsrel`.
///
/// The test returns [`ControlFlow::Break`] if the following condition
/// is achieved,
///
/// $$|x_1 - x_0| < \epsabs + \epsrel \\, |x_1|$$
///
/// and returns [`ControlFlow::Continue`] otherwise.
#[doc(alias = "gsl_root_test_delta")]
pub fn test_delta(x1: f64, x0: f64, epsabs: f64, epsrel: f64) -> ControlFlow<()> {
    Error::control_flow(unsafe { sys::gsl_root_test_delta(x1, x0, epsabs, epsrel) })
}

/// Test the residual value `f` against the absolute error bound
/// `epsabs`.
///
/// The test returns [`ControlFlow::Break`] if the following condition
/// is achieved,
///
/// $$|f| < \epsabs$$
///
/// and returns [`ControlFlow::Continue`] otherwise.  This criterion
/// is suitable for situations where the precise location of the root,
/// $x$, is unimportant provided a value can be found where the
/// residual, $|f(x)|$, is small enough.
#[doc(alias = "gsl_root_test_residual")]
pub fn test_residual(f: f64, epsabs: f64) -> ControlFlow<()> {
    Error::control_flow(unsafe { sys::gsl_root_test_residual(f, epsabs) })
}

#[cfg(any(test, doctest))]
mod test {
    /// This doc block will be used to ensure that the closure can't
    /// be set everywhere!
    ///
    /// ```compile_fail
    /// use crate::rgsl::*;
    ///
    /// fn set(root: &mut RootFSolver) {
    ///     let y = "lalal".to_owned();
    ///     root.set(
    ///         |x| { println!("==> {:?}", y); x * x - 5. },
    ///         0.0, 5.0);
    /// }
    ///
    /// let mut root = RootFSolver::brent();
    /// set(&mut root);
    /// ```
    ///
    /// Same but a working version:
    ///
    /// ```
    /// use crate::rgsl::*;
    ///
    /// fn set(root: &mut RootFSolver) {
    ///     root.set(|x| x * x - 5., 0.0, 5.0);
    /// }
    ///
    /// let mut root = RootFSolver::brent();
    /// set(&mut root);
    /// let status = root.iterate();
    /// ```
    use super::*;
    use crate::roots::{test_delta, test_interval};

    /// Test whether stored pointers are valid long enough.
    #[test]
    fn test_pointer() -> Result<(), Error> {
        #[inline(never)]
        fn create<'a>() -> Result<RootFSolver<'a>, Error> {
            let mut s = RootFSolver::brent();
            s.set(|x| x - 2., 0., 10.)?;
            Ok(s)
        }
        let mut s = create()?;
        s.iterate()?;
        Ok(())
    }

    #[test]
    fn test_names() {
        let s = RootFSolver::bisection();
        assert_eq!(s.name(), Type::Bisection);
        let s = RootFSolver::falsepos();
        assert_eq!(s.name(), Type::FalsePos);
        let s = RootFSolver::brent();
        assert_eq!(s.name(), Type::Brent);
    }

    /// Test whether stored pointers are valid long enough.
    #[test]
    fn test_pointer_fdf() -> Result<(), Error> {
        #[inline(never)]
        fn create<'a>() -> Result<RootFdfSolver<'a>, Error> {
            let mut s = RootFdfSolver::newton();
            let f = |x: f64| x - 2.;
            let df = |_: f64| 1.;
            s.set((f, df), 0.)?;
            Ok(s)
        }
        let mut s = create()?;
        s.iterate()?;
        Ok(())
    }

    #[test]
    fn test_names_fdf() {
        let s = RootFdfSolver::newton();
        assert_eq!(s.name(), TypeFdf::Newton);
        let s = RootFdfSolver::secant();
        assert_eq!(s.name(), TypeFdf::Secant);
        let s = RootFdfSolver::steffensen();
        assert_eq!(s.name(), TypeFdf::Steffensen);
    }

    // support functions
    fn quadratic_test_fn(x: f64) -> f64 {
        x.powi(2) - 5.0
    }

    fn quadratic_test_fn_df(x: f64) -> f64 {
        2.0 * x
    }

    fn quadratic_test_fn_fdf(x: f64) -> (f64, f64) {
        (x.powi(2) - 5.0, 2.0 * x)
    }

    #[test]
    fn test_root() -> Result<(), Error> {
        let mut root = RootFSolver::brent();
        root.set(quadratic_test_fn, 0.0, 5.0)?;

        let max_iter = 10;
        let epsabs = 1e-4;
        let epsrel = 1e-7;

        println!("Testing: {:?}", root.name());

        println!("iter, \t [x_lo, x_hi], \t min, \t error");
        for iter in 0..max_iter {
            root.iterate()?;

            // test for convergence
            let r = root.root();
            let x_lo = root.x_lower();
            let x_hi = root.x_upper();

            println!(
                "{} \t [{:.5}, {:.5}] \t {:.5} \t {:.5}",
                iter,
                x_lo,
                x_hi,
                r,
                x_hi - x_lo
            );
            let status = test_interval(x_lo, x_hi, epsabs, epsrel);

            if status.is_break() {
                println!("Converged");
                break;
            }
        }
        Ok(())
    }

    #[test]
    fn test_root_fdf() -> Result<(), Error> {
        // guess value
        let r_expected = 5.0_f64.sqrt();
        let guess_value = 1.0;

        // setup solver
        let mut root = RootFdfSolver::steffensen();
        root.set(
            (
                quadratic_test_fn,
                quadratic_test_fn_df,
                quadratic_test_fn_fdf,
            ),
            guess_value,
        )?;

        // set up iterations
        let max_iter = 20;
        let epsabs = 1e-4;
        let epsrel = 1e-7;

        println!("Testing: {:?}", root.name());

        println!("iter, \t root, \t rel error \t abs error");

        let mut x = guess_value;
        let mut converged = false;
        for iter in 0..max_iter {
            root.iterate()?;

            // test for convergence
            let x_0 = x;
            x = root.root();
            // print results
            println!(
                "{} \t {:.5} \t {:.5} \t {:.5}",
                iter,
                x,
                x - x_0,
                x - r_expected
            );
            // check if iteration succeeded
            if test_delta(x, x_0, epsabs, epsrel).is_break() {
                println!("Converged");
                converged = true;
                break;
            }
        }
        assert!(converged);
        Ok(())
    }
}
