//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Polynomials

This module functions for evaluating and solving polynomials. There
are routines for finding real and complex roots of quadratic and cubic
equations using analytic methods. An iterative polynomial solver is
also available for finding the roots of general polynomials with real
coefficients (of any order).

## References and Further Reading

The balanced-QR method and its error analysis are described in the following papers,

R.S. Martin, G. Peters and J.H. Wilkinson, “The QR Algorithm for Real Hessenberg Matrices”, Numerische Mathematik, 14 (1970), 219–231.
B.N. Parlett and C. Reinsch, “Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors”, Numerische Mathematik, 13 (1969), 293–304.
A. Edelman and H. Murakami, “Polynomial roots from companion matrix eigenvalues”, Mathematics of Computation, Vol. 64, No. 210 (1995), 763–776.
The formulas for divided differences are given in the following texts,

Abramowitz and Stegun, Handbook of Mathematical Functions, Sections 25.1.4 and 25.2.26.
R. L. Burden and J. D. Faires, Numerical Analysis, 9th edition, ISBN 0-538-73351-9, 2011.
*/

use std::borrow::Cow;

#[cfg(feature = "complex")]
use crate::complex::{FromC, ToC};
use crate::Error;
#[cfg(feature = "complex")]
use num_complex::Complex;

pub struct Poly<'a, T>(&'a [T]);

impl<T> Poly<'_, T> {
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl Poly<'_, f64> {
    /// Return the value of the polynomial $P(x) = c_0 + c_1 x + c_2
    /// x^2 + \dots + c_{n-1} x^{n-1}$ where $n$ = `self.len()`, using
    /// Horner’s method.
    #[doc(alias = "gsl_poly_eval")]
    pub fn eval(&self, x: f64) -> f64 {
        unsafe { sys::gsl_poly_eval(self.0.as_ptr(), self.len() as i32, x) }
    }

    /// Return the value of the polynomial with real coefficients
    /// for the complex variable `z`.
    #[doc(alias = "gsl_poly_complex_eval")]
    #[cfg(feature = "complex")]
    pub fn complex_eval(&self, z: &Complex<f64>) -> Complex<f64> {
        unsafe { sys::gsl_poly_complex_eval(self.0.as_ptr(), self.len() as i32, z.unwrap()).wrap() }
    }

    /// This function evaluates a polynomial and its derivatives
    /// storing the results in the array `res`.  The output array
    /// contains the values of $d^k P/d x^k$ for the specified value
    /// of `x` starting with $k = 0$.
    #[doc(alias = "gsl_poly_eval_derivs")]
    pub fn eval_derivs(c: &[f64], x: f64, res: &mut [f64]) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_poly_eval_derivs(
                c.as_ptr(),
                c.len() as _,
                x,
                res.as_mut_ptr(),
                res.len() as _,
            )
        };
        Error::handle(ret, ())
    }
}

#[cfg(feature = "complex")]
impl Poly<'_, Complex<f64>> {
    /// This function evaluates a polynomial with complex coefficients
    /// for the complex variable z.
    #[doc(alias = "gsl_complex_poly_complex_eval")]
    pub fn eval(&self, z: &Complex<f64>) -> Complex<f64> {
        let c: *const Complex<f64> = self.0.as_ptr();
        unsafe {
            sys::gsl_complex_poly_complex_eval(c as *const _, self.len() as i32, z.unwrap()).wrap()
        }
    }
}

/// The functions described here manipulate polynomials stored in Newton’s divided-difference representation. The use of divided-differences
/// is described in Abramowitz & Stegun sections 25.1.4 and 25.2.26, and Burden and Faires, chapter 3, and discussed briefly below.
///
/// Given a function f(x), an nth degree interpolating polynomial P_{n}(x) can be constructed which agrees with f at n+1 distinct points x_0,
/// x_1,...,x_{n}. This polynomial can be written in a form known as Newton’s divided-difference representation:
///
/// P_n(x) = f(x_0) + \sum_(k=1)^n [x_0,x_1,...,x_k] (x-x_0)(x-x_1)...(x-x_(k-1))
///
/// where the divided differences [x_0,x_1,...,x_k] are defined in section 25.1.4 of Abramowitz and Stegun. Additionally, it is possible to
/// construct an interpolating polynomial of degree 2n+1 which also matches the first derivatives of f at the points x_0,x_1,...,x_n. This is
/// called the Hermite interpolating polynomial and is defined as
///
/// H_(2n+1)(x) = f(z_0) + \sum_(k=1)^(2n+1) [z_0,z_1,...,z_k] (x-z_0)(x-z_1)...(x-z_(k-1))
///
/// where the elements of z = \{x_0,x_0,x_1,x_1,...,x_n,x_n\} are defined by z_{2k} = z_{2k+1} = x_k. The divided-differences [z_0,z_1,...,z_k]
/// are discussed in Burden and Faires, section 3.4.
pub struct PolyDD<'a> {
    dd: Cow<'a, [f64]>,
    x: Cow<'a, [f64]>,
}

impl<'a> PolyDD<'a> {
    pub fn len(&self) -> usize {
        self.dd.len()
    }

    pub fn is_empty(&self) -> bool {
        self.dd.is_empty()
    }

    pub fn new(x: &[f64], y: &[f64]) -> Result<PolyDD<'static>, Error> {
        if x.len() != y.len() {
            panic!(
                "rgsl::polynomial::PolyDD::new: \
                x.len() = {} ≠ y.len() = {}",
                x.len(),
                y.len()
            );
        }
        let mut dd = x.to_owned();
        let ret = unsafe {
            sys::gsl_poly_dd_init(dd.as_mut_ptr(), x.as_ptr(), y.as_ptr(), dd.len() as _)
        };
        let dd = PolyDD {
            dd: Cow::Owned(dd),
            x: Cow::Owned(x.to_owned()),
        };
        Error::handle(ret, dd)
    }

    /// This function computes a divided-difference representation of
    /// the interpolating polynomial for the points $(x_i, y_i)$
    /// stored in the arrays `x` and `ya`.  On output the
    /// divided-differences of $(x_i, y_i)$ are stored in the array
    /// `dd`, of the same length as `x` and `y`.  Using the notation
    /// above, `dd[k] = [x_0,x_1,...,x_k]`.
    #[doc(alias = "gsl_poly_dd_init")]
    pub fn init(dd: &'a mut [f64], x: &'a [f64], y: &[f64]) -> Result<Self, Error> {
        assert_eq!(dd.len(), x.len());
        assert_eq!(x.len(), y.len());
        let ret = unsafe {
            sys::gsl_poly_dd_init(dd.as_mut_ptr(), x.as_ptr(), y.as_ptr(), dd.len() as _)
        };
        let dd = PolyDD {
            dd: Cow::Borrowed(dd),
            x: Cow::Borrowed(x),
        };
        Error::handle(ret, dd)
    }

    /// This function evaluates the polynomial stored
    /// divided-difference form in `self` at the point `x`.
    #[doc(alias = "gsl_poly_dd_eval")]
    pub fn eval(&self, x: f64) -> f64 {
        unsafe { sys::gsl_poly_dd_eval(self.dd.as_ptr(), self.x.as_ptr(), self.len() as _, x) }
    }

    /// This function converts the divided-difference representation
    /// of a polynomial to a Taylor expansion.  On output the Taylor
    /// coefficients of the polynomial expanded about the point `xp`
    /// are stored in the array `c` of length `self.len()`.  A
    /// workspace of the same length must be provided in the array
    /// `w`.
    #[doc(alias = "gsl_poly_dd_taylor")]
    pub fn taylor(
        &self,
        c: &mut [f64],
        xp: f64,
        w: &mut [f64], // TODO: make it optional
    ) -> Result<(), Error> {
        assert_eq!(self.len(), c.len());
        assert_eq!(self.len(), w.len());
        let ret = unsafe {
            sys::gsl_poly_dd_taylor(
                c.as_mut_ptr(),
                xp,
                self.dd.as_ptr(),
                self.x.as_ptr(),
                self.len() as _,
                w.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }

    /// This function computes a divided-difference representation of
    /// the interpolating Hermite polynomial for the points $(x_i,
    /// y_i)$ stored in the arrays `x` and `y`.  Hermite interpolation
    /// constructs polynomials which also match first derivatives
    /// $dy/dx$ which are provided in the array `dy` of the same
    /// length.  The first derivatives can be incorported into the
    /// usual divided-difference algorithm by forming a new dataset $z
    /// = \{x_0,x_0,x_1,x_1,...\}$, which is stored in the array `z`
    /// of length `2 * x.len()` on output.  On output the
    /// divided-differences of the Hermite representation are stored
    /// in the array `dd`, also of length `2 * x.len()`.  Using the
    /// notation above, `dd[k] = [z_0,z_1,...,z_k]`.  The resulting
    /// Hermite polynomial can be evaluated by calling
    /// [`PolyDD::eval`] on the returned [`PolyDD`].
    #[doc(alias = "gsl_poly_dd_hermite_init")]
    pub fn poly_dd_hermite_init(
        dd: &'a mut [f64],
        z: &'a mut [f64],
        x: &[f64],
        y: &[f64],
        dy: &[f64],
    ) -> Result<Self, Error> {
        assert_eq!(x.len(), y.len());
        assert_eq!(y.len(), dy.len());
        let ret = unsafe {
            sys::gsl_poly_dd_hermite_init(
                dd.as_mut_ptr(),
                z.as_mut_ptr(),
                x.as_ptr(),
                y.as_ptr(),
                dy.as_ptr(),
                dd.len() as _,
            )
        };
        let dd = PolyDD {
            dd: Cow::Borrowed(dd),
            x: Cow::Borrowed(z),
        };
        Error::handle(ret, dd)
    }
}

/// Represent $a x^2 + b x + c$.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quadratic {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuadraticRoots<T> {
    None,
    One(T),
    Two(T, T),
}

impl Quadratic {
    /// Return the real roots of the quadratic equation,
    ///
    /// $$a x^2 + b x + c = 0.$$
    ///
    /// When two real roots are found they are returned in ascending
    /// order.  A single root is returned if `a` = 0.  The case of
    /// coincident roots is not considered special.  For example
    /// $(x-1)^2=0$ will have two roots, which happen to have exactly
    /// equal values.
    ///
    /// The number of roots found depends on the sign of the discriminant
    /// $b^2 - 4 a c$.  This will be subject to rounding and cancellation
    /// errors when computed in double precision, and will also be subject
    /// to errors if the coefficients of the polynomial are inexact. These
    /// errors may cause a discrete change in the number of
    /// roots. However, for polynomials with small integer coefficients
    /// the discriminant can always be computed exactly.
    ///
    /// # Example
    ///
    /// ```
    /// use rgsl::polynomials::{Quadratic, QuadraticRoots};
    /// let r = Quadratic { a: 1., b: -2., c: 1. }.real_roots();
    /// assert_eq!(r, QuadraticRoots::Two(1., 1.));
    /// let r = Quadratic { a: 0., b: 0., c: 0. }.real_roots();
    /// assert_eq!(r, QuadraticRoots::None);
    /// ```
    #[doc(alias = "gsl_poly_solve_quadratic")]
    pub fn real_roots(&self) -> QuadraticRoots<f64> {
        let mut x0 = 0.;
        let mut x1 = 0.;
        let n = unsafe { sys::gsl_poly_solve_quadratic(self.a, self.b, self.c, &mut x0, &mut x1) };
        match n {
            0 => QuadraticRoots::None,
            1 => QuadraticRoots::One(x0),
            2 => QuadraticRoots::Two(x0, x1),
            _ => unreachable!(),
        }
    }

    /// Return the complex roots of the quadratic equation
    /// $a z^2 + b z + c = 0$.
    ///
    /// The roots are returned in ascending order, sorted first by
    /// their real components and then by their imaginary components.
    ///
    /// # Example
    ///
    /// ```
    /// use rgsl::polynomials::{Quadratic, QuadraticRoots};
    /// use num_complex::Complex;
    /// let r = Quadratic { a: 1., b: -2., c: 1. }.roots();
    /// let one = Complex::new(1., 0.);
    /// assert_eq!(r, QuadraticRoots::Two(one, one));
    /// let r = Quadratic { a: 0., b: 0., c: 0. }.roots();
    /// assert_eq!(r, QuadraticRoots::None);
    /// ```
    #[doc(alias = "gsl_poly_complex_solve_quadratic")]
    #[cfg(feature = "complex")]
    pub fn roots(&self) -> QuadraticRoots<Complex<f64>> {
        let mut z0 = Complex::new(0., 0.);
        let mut z1 = Complex::new(0., 0.);
        let n = unsafe {
            sys::gsl_poly_complex_solve_quadratic(
                self.a,
                self.b,
                self.c,
                (&mut z0).unwrap(),
                (&mut z1).unwrap(),
            )
        };
        match n {
            0 => QuadraticRoots::None,
            1 => QuadraticRoots::One(z0),
            2 => QuadraticRoots::Two(z0, z1),
            _ => unreachable!(),
        }
    }
}

/// Represent $x^3 + a x^2 + b x + c$.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cubic {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CubicRoots<T> {
    One(T),
    Three(T, T, T),
}

impl Cubic {
    /// Return the real roots of the cubic equation,
    /// $x^3 + a x^2 + b x + c = 0$.
    ///
    /// with a leading coefficient of unity.  The real roots returned
    /// in ascending order.  The case of coincident roots is not
    /// considered special.  For example, the equation $(x-1)^3=0$
    /// will have three roots with exactly equal values.  As in the
    /// quadratic case, finite precision may cause equal or
    /// closely-spaced real roots to move off the real axis into the
    /// complex plane, leading to a discrete change in the number of
    /// real roots.
    #[doc(alias = "gsl_poly_solve_cubic")]
    pub fn real_roots(&self) -> CubicRoots<f64> {
        let mut x0 = 0.;
        let mut x1 = 0.;
        let mut x2 = 0.;
        let n =
            unsafe { sys::gsl_poly_solve_cubic(self.a, self.b, self.c, &mut x0, &mut x1, &mut x2) };
        match n {
            1 => CubicRoots::One(x0),
            3 => CubicRoots::Three(x0, x1, x2),
            _ => unreachable!(),
        }
    }

    /// Return the complex roots of the cubic equation
    /// ^z^3 + a z^2 + b z + c = 0$.
    ///
    /// The complex roots are returned in ascending order, sorted
    /// first by their real components and then by their imaginary
    /// components.
    #[doc(alias = "gsl_poly_complex_solve_cubic")]
    #[allow(unknown_lints, clippy::useless_transmute)]
    #[cfg(feature = "complex")]
    pub fn roots(&self) -> (Complex<f64>, Complex<f64>, Complex<f64>) {
        let mut z0 = Complex::new(0., 0.);
        let mut z1 = Complex::new(0., 0.);
        let mut z2 = Complex::new(0., 0.);
        let _n = unsafe {
            sys::gsl_poly_complex_solve_cubic(
                self.a,
                self.b,
                self.c,
                (&mut z0).unwrap(),
                (&mut z1).unwrap(),
                (&mut z2).unwrap(),
            )
        };
        (z0, z1, z2)
    }
}
