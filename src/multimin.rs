//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Multidimensional Minimization

This module contains routines for finding minima of arbitrary
multidimensional functions.  It provides low level components for a
variety of iterative minimizers and convergence tests.  These can be
combined by the user to achieve the desired solution, while providing
full access to the intermediate steps of the algorithms.  Each class
of methods uses the same framework, so that you can switch between
minimizers at runtime without needing to recompile your program.  Each
instance of a minimizer keeps track of its own state, allowing the
minimizers to be used in multi-threaded programs.  The minimization
algorithms can be used to maximize a function by inverting its sign.

## Overview

The problem of multidimensional minimization requires finding a point
$x$ such that the scalar function,

$$f(x_1,..., x_n)$$

takes a value which is lower than at any neighboring point.  For
smooth functions the gradient $g = \nabla f$ vanishes at the minimum.
In general there are no bracketing methods available for the
minimization of n-dimensional functions.  The algorithms proceed from
an initial guess using a search algorithm which attempts to move in a
downhill direction.

Algorithms making use of the gradient of the function perform a
one-dimensional line minimisation along this direction until the
lowest point is found to a suitable tolerance.  The search direction
is then updated with local information from the function and its
derivatives, and the whole process repeated until the true
n-dimensional minimum is found.

Algorithms which do not require the gradient of the function use
different strategies.  For example, the Nelder-Mead Simplex algorithm
maintains $n+1$ trial parameter vectors as the vertices of a
$n$-dimensional simplex.  On each iteration it tries to improve the
worst vertex of the simplex by geometrical transformations.  The
iterations are continued until the overall size of the simplex has
decreased sufficiently.

Both types of algorithms use a standard framework.  The user provides
a high-level driver for the algorithms, and the library provides the
individual functions necessary for each of the steps.  There are three
main phases of the iteration. The steps are,

- initialize minimizer state, `s`, for algorithm `t`,
- update `s` using the iteration `t`,
- test `s` for convergence, and repeat iteration if necessary.

Each iteration step consists either of an improvement to the
line-minimisation in the current direction or an update to the search
direction itself.  The state for the minimizers is held in
[`MinimizerFdf`] or [`Minimizer`].

## Caveats

Note that the minimization algorithms can only search for one local
minimum at a time.  When there are several local minima in the search
area, the first minimum to be found will be returned; however it is
difficult to predict which of the minima this will be.  In most cases,
no error will be reported if you try to find a local minimum in an
area where there is more than one.

It is also important to note that the minimization algorithms find
local minima; there is no way to determine whether a minimum is a
global minimum of the function in question.

## Example

TODO
 */

use crate::{
    Error,
    ffi::FFI,
    vector::{AsVector, VecF64},
    view::{AsView, View},
};
use std::{ffi::c_void, ops::ControlFlow};

/// Type of minimizer without derivatives.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Type {
    Simplex,
    Simplex2,
    Simplex2Rand,
}

impl Type {
    #[inline]
    fn to_c(self) -> *const sys::gsl_multimin_fminimizer_type {
        unsafe {
            match self {
                Type::Simplex => sys::gsl_multimin_fminimizer_nmsimplex,
                Type::Simplex2 => sys::gsl_multimin_fminimizer_nmsimplex2,
                Type::Simplex2Rand => sys::gsl_multimin_fminimizer_nmsimplex2rand,
            }
        }
    }
}

ffi_wrapper!(
    /// Minimization without derivatives.
    ///
    /// # Examples
    ///
    /// Minimize the function $(x₀, x₁) ↦ (x₀ - 1)² + (x₁ - 1)²$ using
    /// [`VecF64`] to represent vectors.
    ///
    /// ```
    /// use rgsl::{VecF64, multimin::Minimizer};
    /// let x = VecF64::from_slice(&[0., 0.]);
    /// let step_size = VecF64::from_slice(&[0.1, 0.1]);
    /// let mut solver = Minimizer::simplex(x.len());
    /// let f = |x: &VecF64| (x.get(0) - 1.).powi(2) + (x.get(1) - 1.).powi(2);
    /// solver.set(f, &x, &step_size)?;
    /// for _ in 0 .. 120 {
    ///     solver.iterate()?;
    /// }
    /// assert_eq!(solver.x(), VecF64::from_slice(&[1., 1.]));
    /// # Ok::<(), rgsl::Error>(())
    /// ```
    ///
    /// Here is the minimization of the same function using slices
    /// `&[f64]` to represent vectors.
    ///
    /// ```
    /// use rgsl::{VecF64, multimin::Minimizer};
    /// let x = [0., 0.].as_slice();
    /// let step_size = [0.1, 0.1].as_slice();
    /// let mut solver = Minimizer::simplex(x.len());
    /// let f = |x: &[f64]| (x[0] - 1.).powi(2) + (x[1] - 1.).powi(2);
    /// solver.set(f, x, step_size)?;
    /// for _ in 0 .. 120 {
    ///     solver.iterate()?;
    /// }
    /// assert_eq!(solver.x(), &[1., 1.]);
    /// # Ok::<(), rgsl::Error>(())
    /// ```
    ///
    /// With the feature `ndarray`, `Array1<f64>` may also be used for
    /// vectors.  To enable it for you vector type, implement the
    /// trait [`AsVector`].
    Minimizer<'a, V: AsVector + ?Sized>,
    *mut sys::gsl_multimin_fminimizer,
    gsl_multimin_fminimizer_free
    // These fields hold the GSL structure and function with stable address.
    ;f_struct: Option<Box<sys::gsl_multimin_function_struct>> => None;
    ;f: Option<Box<dyn FnMut(&V) -> f64 + 'a>> => None;
);

impl<'a, V: AsVector + ?Sized> Minimizer<'a, V> {
    /// Creates a minimizer of type `t` for an `n`-dimensional function.
    /// Panic if there is insufficient memory.
    #[doc(alias = "gsl_multimin_fminimizer_alloc")]
    pub fn new(t: Type, n: usize) -> Self {
        let ptr = unsafe { sys::gsl_multimin_fminimizer_alloc(t.to_c(), n) };
        if ptr.is_null() {
            panic!("rgsl::multimin::Minimizer::new: out of memory");
        }
        Self::wrap(ptr)
    }

    pub fn simplex(n: usize) -> Self {
        Self::new(Type::Simplex, n)
    }

    pub fn simplex2(n: usize) -> Self {
        Self::new(Type::Simplex2, n)
    }

    pub fn simplex2_rand(n: usize) -> Self {
        Self::new(Type::Simplex2Rand, n)
    }

    /// Initialize the minimizer to minimize the function `f`,
    /// starting from the initial point `x`.  The size of the initial
    /// trial steps is given in vector `step_size`. The precise
    /// meaning of this parameter depends on the method used.
    #[doc(alias = "gsl_multimin_fminimizer_set")]
    pub fn set<F: FnMut(&V) -> f64 + 'a>(
        &mut self,
        f: F,
        x: &V,
        step_size: &V,
    ) -> Result<(), Error> {
        unsafe extern "C" fn inner_f<V: AsVector + ?Sized, F: FnMut(&V) -> f64>(
            x: *const sys::gsl_vector,
            params: *mut c_void,
        ) -> f64 {
            let f = unsafe { &mut *params.cast::<F>() };
            let vx = unsafe { V::view_from_ptr(x) };
            f(&vx)
        }
        let mut f = Box::new(f);
        let mut f_struct = Box::new(sys::gsl_multimin_function_struct {
            f: Some(inner_f::<V, F>),
            n: V::len(x),
            params: &mut *f as *mut F as *mut _,
        });
        self.f = Some(f);

        let x = V::as_gsl_vector(x);
        let step_size = V::as_gsl_vector(step_size);
        let ret = unsafe {
            sys::gsl_multimin_fminimizer_set(
                self.unwrap_unique(),
                &mut *f_struct,
                &*x,         // Copied inside the GSL value in `self`
                &*step_size, // Only used by this C fn
            )
        };
        self.f_struct = Some(f_struct);
        Error::handle(ret, ())
    }

    /// This function returns the type of the minimizer.
    #[doc(alias = "gsl_multimin_fminimizer_name")]
    pub fn name(&self) -> Type {
        let n = unsafe { sys::gsl_multimin_fminimizer_name(self.unwrap_shared()) };
        map_name!(
            rgsl::multimin::Minimizer,
            [
                (c"nmsimplex", Type::Simplex),
                (c"nmsimplex2", Type::Simplex2),
                (c"nmsimplex2rand", Type::Simplex2Rand)
            ],
            n,
            Type
        )
    }

    /// Location of the minimum at the current point of iteration.
    #[doc(alias = "gsl_multimin_fminimizer_x")]
    pub fn x(&self) -> V::View<'_> {
        unsafe {
            let ptr = sys::gsl_multimin_fminimizer_x(self.unwrap_shared());
            V::view_from_ptr(ptr)
        }
    }

    /// Value of the minimum at the current point of iteration.
    #[doc(alias = "gsl_multimin_fminimizer_minimum")]
    pub fn minimum(&self) -> f64 {
        unsafe { sys::gsl_multimin_fminimizer_minimum(self.unwrap_shared()) }
    }

    /// Return the specific characteristic size for the minimizer.
    #[doc(alias = "gsl_multimin_fminimizer_size")]
    pub fn size(&self) -> f64 {
        unsafe { sys::gsl_multimin_fminimizer_size(self.unwrap_shared()) }
    }

    /// Perform a single iteration of the minimizer `self`.  If the
    /// iteration encounters an unexpected problem then an error code
    /// will be returned.  The error code [`Error::NoProgress`]
    /// signifies that the minimizer is unable to improve on its
    /// current estimate, either due to numerical difficulty or
    /// because a genuine local minimum has been reached.
    #[doc(alias = "gsl_multimin_fminimizer_iterate")]
    pub fn iterate(&mut self) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_multimin_fminimizer_iterate(self.unwrap_unique()) };
        Error::handle(ret, ())
    }
}

/// Specify which data-types may be used in minimization with derivatives.
pub trait Fdf<V: AsVector + ?Sized> {
    /// Return the value of the function $f$ to be minimized.
    fn f(&mut self, x: &V) -> f64;
    /// Set `g` to the gradient of $f$: `g[i]`$= ∂f/∂xᵢ$.
    fn df(&mut self, x: &V, g: &mut V);
    /// Return the value of the function $f$ to be minimized and set
    /// `g` to its gradient.
    fn fdf(&mut self, x: &V, g: &mut V) -> f64;
}

impl<V, F, G> Fdf<V> for (F, G)
where
    V: AsVector + ?Sized,
    F: FnMut(&V) -> f64,
    G: FnMut(&V, &mut V),
{
    #[inline]
    fn f(&mut self, x: &V) -> f64 {
        self.0(x)
    }
    #[inline]
    fn df(&mut self, x: &V, g: &mut V) {
        self.1(x, g)
    }
    #[inline]
    fn fdf(&mut self, x: &V, g: &mut V) -> f64 {
        self.1(x, g);
        self.0(x)
    }
}

impl<V, FG> Fdf<V> for FG
where
    V: AsVector + ?Sized,
    FG: FnMut(&V, Option<&mut V>) -> f64,
{
    #[inline]
    fn f(&mut self, x: &V) -> f64 {
        self(x, None)
    }
    #[inline]
    fn df(&mut self, x: &V, g: &mut V) {
        self(x, Some(g));
    }
    #[inline]
    fn fdf(&mut self, x: &V, g: &mut V) -> f64 {
        self(x, Some(g))
    }
}

impl<V, F, G, FG> Fdf<V> for (F, G, FG)
where
    V: AsVector + ?Sized,
    F: FnMut(&V) -> f64,
    G: FnMut(&V, &mut V),
    FG: FnMut(&V, &mut V) -> f64,
{
    #[inline]
    fn f(&mut self, x: &V) -> f64 {
        self.0(x)
    }
    #[inline]
    fn df(&mut self, x: &V, g: &mut V) {
        self.1(x, g)
    }
    #[inline]
    fn fdf(&mut self, x: &V, g: &mut V) -> f64 {
        self.2(x, g)
    }
}

/// Minimization algorithm using gradients.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TypeFdf {
    /// Fletcher-Reeves conjugate gradient algorithm.  See
    /// [`MinimizerFdf::conjugate_fr`] for more details.
    Conjugate_fr,
    /// Polak-Ribiere conjugate gradient algorithm.  See
    /// [`MinimizerFdf::conjugate_pr`] for more details.
    Conjugate_pr,
    /// Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm.  See
    /// [`MinimizerFdf::bfgs`] for more details.
    BFGS,
    /// Most efficient version of the Broyden-Fletcher-Goldfarb-Shanno
    /// (BFGS) algorithm.  See [`MinimizerFdf::bfgs2`] for more details.
    BFGS2,
    /// Steepest descent algorithm.  See
    /// [`MinimizerFdf::steepest_descent`] for more details.
    Steepest_descent,
}

impl TypeFdf {
    #[inline]
    fn to_c(self) -> *const sys::gsl_multimin_fdfminimizer_type {
        unsafe {
            match self {
                TypeFdf::Conjugate_fr => sys::gsl_multimin_fdfminimizer_conjugate_fr,
                TypeFdf::Conjugate_pr => sys::gsl_multimin_fdfminimizer_conjugate_pr,
                TypeFdf::BFGS => sys::gsl_multimin_fdfminimizer_vector_bfgs,
                TypeFdf::BFGS2 => sys::gsl_multimin_fdfminimizer_vector_bfgs2,
                TypeFdf::Steepest_descent => sys::gsl_multimin_fdfminimizer_steepest_descent,
            }
        }
    }
}

ffi_wrapper!(
    /// Minimization using derivatives.
    ///
    /// # Example
    ///
    /// Minimize the function $(x₀, x₁) ↦ (x₀ - 1)² + (x₁ - 2)²$ using
    /// [`VecF64`] to represent vectors.
    ///
    /// ```
    /// use rgsl::{Error, VecF64, multimin::{self as mmin, MinimizerFdf}};
    /// use std::ops::ControlFlow;
    /// let x = VecF64::from_slice(&[0., 0.]);
    /// let mut solver = MinimizerFdf::bfgs2(x.len());
    /// let f = |x: &VecF64| (x.get(0) - 1.).powi(2) + (x.get(1) - 2.).powi(2);
    /// let df = |x: &VecF64, g: &mut VecF64| {
    ///     g.set(0, 2. * (x.get(0) - 1.));
    ///     g.set(1, 2. * (x.get(1) - 2.));
    /// };
    /// solver.set((f, df), &x, 0.1, 1e-3)?;
    /// for _ in 0 .. 100 {
    ///     solver.iterate()?;
    ///     let status = mmin::test_gradient(&*solver.gradient(), 1e-6);
    ///     if matches!(status, ControlFlow::Break(())) {
    ///         break
    ///     }
    /// }
    /// assert!((solver.x().get(0) - 1.).abs() < 1e-12);
    /// assert!((solver.x().get(1) - 2.).abs() < 1e-12);
    /// # Ok::<(), rgsl::Error>(())
    /// ```
    MinimizerFdf<'a, V: AsVector + ?Sized>,
    *mut sys::gsl_multimin_fdfminimizer,
    gsl_multimin_fdfminimizer_free
    ;fdf_struct: Option<Box<sys::gsl_multimin_function_fdf_struct>> => None;
    ;fdf: Option<Box<dyn Fdf<V> + 'a>> => None;
);

impl<'a, V: AsVector + ?Sized> MinimizerFdf<'a, V> {
    /// Creates a minimizer of type `t` for an `n`-dimensional function.
    /// If there is insufficient memory to create the minimizer then
    /// the function returns a `None`.
    #[doc(alias = "gsl_multimin_fdfminimizer_alloc")]
    pub fn new(t: TypeFdf, n: usize) -> Self {
        let ptr = unsafe { sys::gsl_multimin_fdfminimizer_alloc(t.to_c(), n) };
        if ptr.is_null() {
            panic!("rgsl::multimin::MinimizerFdf::new: cannot allocate");
        }
        Self::wrap(ptr)
    }

    /// This is the Fletcher-Reeves conjugate gradient algorithm. The
    /// conjugate gradient algorithm proceeds as a succession of line
    /// minimizations. The sequence of search directions is used to
    /// build up an approximation to the curvature of the function in
    /// the neighborhood of the minimum.
    ///
    /// An initial search direction $p$ is chosen using the gradient,
    /// and line minimization is carried out in that direction.  The
    /// accuracy of the line minimization is specified by the
    /// parameter `tol`.  The minimum along this line occurs when the
    /// function gradient $g$ and the search direction $p$ are
    /// orthogonal.  The line minimization terminates when $p·g < \tol
    /// |p| |g|$.  The search direction is updated using the
    /// Fletcher-Reeves formula $p′ = g′ - β p$ where $β = -|g′|²/|g|²$,
    /// and the line minimization is then repeated for the new search
    /// direction.
    #[doc(alias = "gsl_multimin_fdfminimizer_conjugate_fr")]
    pub fn conjugate_fr(n: usize) -> Self {
        Self::new(TypeFdf::Conjugate_fr, n)
    }

    /// This is the Polak-Ribiere conjugate gradient algorithm.  It is
    /// similar to the Fletcher-Reeves method (see
    /// [`Self::conjugate_fr`]), differing only in the choice of the
    /// coefficient $β$.  Both methods work well when the evaluation
    /// point is close enough to the minimum of the objective function
    /// that it is well approximated by a quadratic hypersurface.
    #[doc(alias = "gsl_multimin_fdfminimizer_conjugate_pr")]
    pub fn conjugate_pr(n: usize) -> Self {
        Self::new(TypeFdf::Conjugate_pr, n)
    }

    /// This method uses the vector Broyden-Fletcher-Goldfarb-Shanno
    /// (BFGS) algorithm.  This is a quasi-Newton method which builds
    /// up an approximation to the second derivatives of the function
    /// $f$ using the difference between successive gradient vectors.
    /// By combining the first and second derivatives the algorithm is
    /// able to take Newton-type steps towards the function minimum,
    /// assuming quadratic behavior in that region.
    #[doc(alias = "gsl_multimin_fdfminimizer_vector_bfgs")]
    pub fn bfgs(n: usize) -> Self {
        Self::new(TypeFdf::BFGS, n)
    }

    /// Version of the BFGS minimizer that is the most efficient
    /// version available, and is a faithful implementation of the
    /// line minimization scheme described in Fletcher’s Practical
    /// Methods of Optimization, Algorithms 2.6.2 and 2.6.4.  It
    /// supersedes the original bfgs routine and requires
    /// substantially fewer function and gradient evaluations.  The
    /// user-supplied tolerance `tol` (see [`Self::set`]) corresponds
    /// to the parameter $σ$ used by Fletcher.  A value of 0.1 is
    /// recommended for typical use (larger values correspond to less
    /// accurate line searches).
    #[doc(alias = "gsl_multimin_fdfminimizer_vector_bfgs2")]
    pub fn bfgs2(n: usize) -> Self {
        Self::new(TypeFdf::BFGS2, n)
    }

    /// The steepest descent algorithm follows the downhill gradient
    /// of the function at each step.  When a downhill step is
    /// successful the step-size is increased by a factor of two.  If
    /// the downhill step leads to a higher function value then the
    /// algorithm backtracks and the step size is decreased using the
    /// parameter `tol`.  A suitable value of tol for most
    /// applications is 0.1.  The steepest descent method is
    /// inefficient and is included only for demonstration purposes.
    #[doc(alias = "gsl_multimin_fdfminimizer_steepest_descent")]
    pub fn steepest_descent(n: usize) -> Self {
        Self::new(TypeFdf::Steepest_descent, n)
    }

    /// Initialize the minimizer to minimize the function (and its
    /// derivative) `fdf` starting from the initial point `x`.
    ///
    /// The argument `fdf` may be:
    /// - a couple of functions `(f, df)` where `f(x: &V) -> f64`
    ///   computes the function to minimize and `df(x: &V, g: &mut g)`
    ///   stores the gradient of `f` at `x` in `g`;
    /// - a function `df(x: &V, g: Option<&mut V>) -> f64` that
    ///   returns the value at `x` of the function to minimize and
    ///   stores its gradient in `g` if provided;
    /// - a triplet `(f, df, fdf)` where `(f, df)` is as in the first
    ///   case and `fdf(x: &V, g: &mut V) -> f64` is as in the second
    ///   point where the gradient storage `g` is provided.
    ///
    /// The size of the first trial step is given by `step_size`.  The
    /// accuracy of the line minimization is specified by `tol`.  The
    /// precise meaning of this parameter depends on the method used.
    /// Typically the line minimization is considered successful if
    /// the gradient of the function $g$ is orthogonal to the current
    /// search direction $p$ to a relative accuracy of `tol`, where $p
    /// · g < \tol |p| |g|$.  A `tol` value of 0.1 is suitable for most
    /// purposes, since line minimization only needs to be carried out
    /// approximately.  Note that setting `tol` to zero will force the
    /// use of “exact” line-searches, which are extremely expensive.
    #[doc(alias = "gsl_multimin_fdfminimizer_set")]
    pub fn set<F: Fdf<V> + 'a>(
        &mut self,
        fdf: F,
        x: &V,
        step_size: f64,
        tol: f64,
    ) -> Result<(), Error> {
        unsafe extern "C" fn inner_f<V: AsVector + ?Sized, F: Fdf<V>>(
            x: *const sys::gsl_vector,
            params: *mut c_void,
        ) -> f64 {
            let t = unsafe { &mut *params.cast::<F>() };
            let vx = unsafe { V::view_from_ptr(x) };
            t.f(&vx)
        }

        unsafe extern "C" fn inner_df<V: AsVector + ?Sized, F: Fdf<V>>(
            x: *const sys::gsl_vector,
            params: *mut c_void,
            g: *mut sys::gsl_vector,
        ) {
            unsafe {
                let t = &mut *params.cast::<F>();
                let vx = V::view_from_ptr(x);
                let mut vg = V::view_from_mut_ptr(g);
                t.df(&vx, &mut vg);
            }
        }
        unsafe extern "C" fn inner_fdf<V: AsVector + ?Sized, F: Fdf<V>>(
            x: *const sys::gsl_vector,
            params: *mut c_void,
            f: *mut f64,
            g: *mut sys::gsl_vector,
        ) {
            unsafe {
                let t = &mut *params.cast::<F>();
                let vx = V::view_from_ptr(x);
                let mut vg = V::view_from_mut_ptr(g);
                *f = t.fdf(&vx, &mut vg);
            }
        }
        let mut fdf = Box::new(fdf);
        let mut fdf_struct = Box::new(sys::gsl_multimin_function_fdf_struct {
            f: Some(inner_f::<V, F>),
            df: Some(inner_df::<V, F>),
            fdf: Some(inner_fdf::<V, F>),
            n: V::len(x),
            params: &mut *fdf as *mut F as *mut _,
        });
        self.fdf = Some(fdf);

        let x = V::as_gsl_vector(x);
        let ret = unsafe {
            sys::gsl_multimin_fdfminimizer_set(
                self.unwrap_unique(),
                &mut *fdf_struct,
                &*x, // Copied inside the GSL value in `self`
                step_size,
                tol,
            )
        };
        self.fdf_struct = Some(fdf_struct);
        Error::handle(ret, ())
    }

    /// This function returns the type of the minimizer.
    #[doc(alias = "gsl_multimin_fdfminimizer_name")]
    pub fn name(&self) -> TypeFdf {
        let n = unsafe { sys::gsl_multimin_fdfminimizer_name(self.unwrap_shared()) };
        map_name!(
            rgsl::multimin::MinimizerFdf,
            [
                (c"conjugate_fr", TypeFdf::Conjugate_fr),
                (c"vector_bfgs", TypeFdf::BFGS),
                (c"vector_bfgs2", TypeFdf::BFGS2),
                (c"steepest_descent", TypeFdf::Steepest_descent),
            ],
            n,
            TypeFdf
        )
    }

    /// Returns the current best estimate of the location of the minimum.
    #[doc(alias = "gsl_multimin_fdfminimizer_x")]
    pub fn x(&self) -> View<'_, VecF64> {
        VecF64::as_view(unsafe { sys::gsl_multimin_fdfminimizer_x(self.unwrap_shared()) })
    }

    /// Returns the value of the function at the minimum.
    #[doc(alias = "gsl_multimin_fdfminimizer_minimum")]
    pub fn minimum(&self) -> f64 {
        unsafe { sys::gsl_multimin_fdfminimizer_minimum(self.unwrap_shared()) }
    }

    /// Returns the gradient of the function at the minimum.
    #[doc(alias = "gsl_multimin_fdfminimizer_gradient")]
    pub fn gradient(&self) -> View<'_, VecF64> {
        VecF64::as_view(unsafe { sys::gsl_multimin_fdfminimizer_gradient(self.unwrap_shared()) })
    }

    /// Returns the last step increment of the estimate.
    #[doc(alias = "gsl_multimin_fdfminimizer_dx")]
    pub fn dx(&self) -> View<'_, VecF64> {
        VecF64::as_view(unsafe { sys::gsl_multimin_fdfminimizer_dx(self.unwrap_shared()) })
    }

    /// Perform a single iteration of the minimizer `self`.  If the
    /// iteration encounters an unexpected problem then an error code
    /// will be returned. The error code `Error::NoProgress` signifies
    /// that the minimizer is unable to improve on its current
    /// estimate, either due to numerical difficulty or because a
    /// genuine local minimum has been reached.
    #[doc(alias = "gsl_multimin_fdfminimizer_iterate")]
    pub fn iterate(&mut self) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_multimin_fdfminimizer_iterate(self.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// This function resets the minimizer to use the current point as
    /// a new starting point.
    #[doc(alias = "gsl_multimin_fdfminimizer_restart")]
    pub fn restart(&mut self) -> i32 {
        unsafe { sys::gsl_multimin_fdfminimizer_restart(self.unwrap_unique()) }
    }
}

/// Test the minimizer specific characteristic size (if applicable to
/// the used minimizer) against absolute tolerance `epsabs`.  The test
/// returns `ControlFlow::Break(())` if the size is smaller than
/// tolerance, otherwise `ControlFlow::Continue(())` is returned.
#[doc(alias = "gsl_multimin_test_size")]
pub fn test_size(size: f64, epsabs: f64) -> ControlFlow<()> {
    Error::control_flow(unsafe { sys::gsl_multimin_test_size(size, epsabs) })
}

/// Test the norm of the gradient `g` against the absolute tolerance
/// `epsabs`.  The gradient of a multidimensional function goes to
/// zero at a minimum.  The test returns `ControlFlow::Break(())` if
/// the following condition is achieved, $|g| < $`epsabs` and returns
/// `ControlFlow::Continue(())` otherwise.  A suitable choice of
/// `epsabs` can be made from the desired accuracy in the function for
/// small variations in `x`.  The relationship between these
/// quantities is given by $δf = g \, δx$.
#[doc(alias = "gsl_multimin_test_gradient")]
pub fn test_gradient<V>(g: &V, epsabs: f64) -> ControlFlow<()>
where
    V: AsVector + ?Sized,
{
    let g = V::as_gsl_vector(g);
    Error::control_flow(unsafe { sys::gsl_multimin_test_gradient(&*g, epsabs) })
}

#[cfg(any(test, doctest))]
mod test {
    /// This doc block will be used to ensure that the closure can't
    /// be set everywhere!
    ///
    /// ```compile_fail
    /// use crate::rgsl::{VecF64, multimin::{Minimizer,Type}};
    ///
    /// fn set(m: &mut Minimizer<VecF64>) {
    ///     let dummy = "lalal".to_owned();
    ///     m.set(|x| {
    ///         println!("==> {:?}", dummy);
    ///         x.get(0) + x.get(1)},
    ///         &VecF64::from_slice(&[-10.0, 1.0]),
    ///         &VecF64::from_slice(&[1.0, 1.0]));
    ///
    ///     let mut mint = Minimizer::new(Type::Simplex, 2);
    ///     set(&mut mint);
    ///     let _status = mint.iterate();
    /// }
    /// ```
    ///
    /// Same but a working version:
    ///
    /// ```
    /// use crate::rgsl::{VecF64, multimin::{Minimizer,Type}};
    ///
    /// fn set(m: &mut Minimizer<VecF64>) {
    ///     m.set(|x| {
    ///         x.get(0) + x.get(1)},
    ///         &VecF64::from_slice(&[-10.0, 1.0]),
    ///         &VecF64::from_slice(&[1.0, 1.0]));
    ///
    ///     let mut mint = Minimizer::new(Type::Simplex, 2);
    ///     set(&mut mint);
    ///     let _status = mint.iterate();
    /// }
    /// ```
    use super::*;

    #[cfg(feature = "ndarray")]
    #[test]
    fn test_multimin_ndarray() -> Result<(), Error> {
        use ndarray::prelude::*;
        let x = Array1::from_vec(vec![0., 0.]);
        let step_size = Array1::from_vec(vec![0.1, 0.1]);
        let mut solver = Minimizer::simplex(x.len());
        let f = |x: &ArrayRef1<f64>| (x[0] - 1.).powi(2) + (x[1] - 1.).powi(2);
        solver.set(f, &x, &step_size)?;
        Ok(())
    }

    #[test]
    fn test_minimizer_pointers() -> Result<(), Error> {
        // This function will move the `Minimizer`, testing that its
        // inner pointers stay valid.
        fn create() -> Result<Minimizer<'static, VecF64>, Error> {
            let mut m = Minimizer::simplex(1);
            let x = VecF64::from_slice(&[1.]);
            let step_size = VecF64::from_slice(&[0.1]);
            m.set(|x| x.get(0), &x, &step_size)?;
            Ok(m)
        }
        let mut m = create()?;
        m.iterate()?;
        Ok(())
    }

    #[test]
    fn test_minimizerfdf_pointers() -> Result<(), Error> {
        // This function will move the `Minimizer`, testing that its
        // inner pointers stay valid.
        fn create() -> Result<MinimizerFdf<'static, VecF64>, Error> {
            let mut m = MinimizerFdf::bfgs(1);
            let x = VecF64::from_slice(&[1.]);
            m.set(
                |x: &VecF64, g: Option<&mut VecF64>| {
                    if let Some(g) = g {
                        g.set(0, 1.);
                    }
                    x.get(0)
                },
                &x,
                0.1,
                1e-3,
            )?;
            Ok(m)
        }
        let mut m = create()?;
        m.iterate()?;
        Ok(())
    }

    fn print_f_state(min: &Minimizer<VecF64>, iter: usize) {
        let f = min.minimum();
        let x = min.x();
        println!(
            "iter: {}, f = {:+.2e}, x = [{:+.5}, {:+.5}]",
            iter,
            f,
            x.get(0),
            x.get(1),
        )
    }

    fn print_fdf_state(min: &MinimizerFdf<VecF64>, iter: usize) {
        let f = min.minimum();
        let x = min.x();
        println!(
            "iter: {}, f = {:+.2e}, x = [{:+.5}, {:+.5}]",
            iter,
            f,
            x.get(0),
            x.get(1),
        )
    }

    const CENTER: (f64, f64) = (1.0, 2.0);
    const SCALE: (f64, f64) = (10.0, 20.0);
    const MINIMUM: f64 = 30.0;

    fn paraboloid(v: &VecF64) -> f64 {
        let x = v.get(0);
        let y = v.get(1);
        let result =
            SCALE.0 * (x - CENTER.0).powf(2.0) + SCALE.1 * (y - CENTER.1).powf(2.0) + MINIMUM;
        result
    }

    #[test]
    fn test_multi_min() -> Result<(), Error> {
        let mut min = Minimizer::new(Type::Simplex2, 2);
        assert_eq!(min.name(), Type::Simplex2);
        let guess_value = VecF64::from_slice(&[5.0, 7.0]);
        let step_size = VecF64::from_slice(&[1.0, 1.0]);

        min.set(paraboloid, &guess_value, &step_size)?;

        let max_iter = 100_usize;
        let eps_abs = 0.01;

        for iter in 0..max_iter {
            // iterate for next value
            min.iterate()?;

            // test for convergence
            let size = min.size();

            // print current iteration
            print_f_state(&min, iter);

            let status = test_size(size, eps_abs);
            if matches!(status, ControlFlow::Break(())) {
                println!("Converged");
                break;
            }
        }
        Ok(())
    }

    #[test]
    fn test_multi_fdf_min() -> Result<(), Error> {
        let mut min = MinimizerFdf::new(TypeFdf::Conjugate_fr, 2);
        assert_eq!(min.name(), TypeFdf::Conjugate_fr);
        let guess_value = VecF64::from_slice(&[5.0, 7.0]);
        let step_size = 0.01;
        let tol = 1e-4;

        fn df(v: &VecF64, g: &mut VecF64) {
            let x = v.get(0);
            let y = v.get(1);
            g.set(0, 2.0 * SCALE.0 * (x - CENTER.0));
            g.set(1, 2.0 * SCALE.1 * (y - CENTER.1));
        }

        fn fdf(v: &VecF64, g: &mut VecF64) -> f64 {
            df(v, g);
            paraboloid(v)
        }

        min.set((paraboloid, df, fdf), &guess_value, step_size, tol)?;

        let max_iter = 100_usize;
        let eps_abs = 0.01;

        for iter in 0..max_iter {
            // iterate for next value
            min.iterate()?;

            // print current iteration
            print_fdf_state(&min, iter);

            let status = test_gradient(&*min.gradient(), eps_abs);
            if matches!(status, ControlFlow::Break(())) {
                println!("Converged");
                break;
            }
        }
        Ok(())
    }
}
