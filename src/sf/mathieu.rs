/*!
# Mathieu Functions

The routines described in this section compute the angular and radial
Mathieu functions, and their characteristic values. Mathieu functions
are the solutions of the following two differential equations:

$$\begin{align}
& \frac{d^2y}{dv^2} + (a - 2q\cos 2v)y = 0
\newline
& \frac{d^2f}{du^2} - (a - 2q\cosh 2u)f = 0
\end{align}$$

The angular Mathieu functions $\ce_n(x,q)$, $\se_n(x,q)$ are the even
and odd periodic solutions of the first equation, which is known as
Mathieu’s equation.  These exist only for the discrete sequence of
characteristic values $a=a_n(q)$ (even-periodic) and $a=b_n(q)$
(odd-periodic).

The radial Mathieu functions $\Mc_n^{(j)}(z,q)$, $\Ms_n^{(j)}(z,q)$
are the solutions of the second equation, which is referred to as
Mathieu’s modified equation.  The radial Mathieu functions of the
first, second, third and fourth kind are denoted by the parameter $j$,
which takes the value 1, 2, 3 or 4.

For more information on the Mathieu functions, see Abramowitz and
Stegun, Chapter 20.
!*/

use crate::ffi::FFI;
use crate::{types, Error};
use std::mem::MaybeUninit;

ffi_wrapper!(
    Mathieu,
    *mut sys::gsl_sf_mathieu_workspace,
    gsl_sf_mathieu_free,
    "Workspace to compute array-based routines."
);

impl Mathieu {
    /// This function returns a workspace for the array versions of
    /// the Mathieu routines.  The arguments `n` and `qmax` specify
    /// the maximum order and q-value of Mathieu functions which can
    /// be computed with this workspace.
    #[doc(alias = "gsl_sf_mathieu_alloc")]
    pub fn new(n: usize, qmax: f64) -> Option<Self> {
        let tmp = unsafe { sys::gsl_sf_mathieu_alloc(n, qmax) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// Return the characteristic values $a_n(q)$ of the Mathieu
    /// function $\ce_n(q,x)$.
    #[doc(alias = "gsl_sf_mathieu_a_e")]
    pub fn a(n: i32, q: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_mathieu_a_e(n, q, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    #[deprecated(since = "8.0.0", note = "Use rgsl::sf::mathieu::Mathieu::a")]
    pub fn mathieu_a(n: i32, q: f64) -> Result<types::Result, Error> {
        Self::a(n, q)
    }

    /// Return the characteristic values $b_n(q)$ of the Mathieu
    /// function $\se_n(q,x)$.
    #[doc(alias = "gsl_sf_mathieu_b_e")]
    pub fn b(n: i32, q: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_mathieu_b_e(n, q, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    #[deprecated(since = "8.0.0", note = "Use rgsl::sf::mathieu::Mathieu::b")]
    pub fn mathieu_b(n: i32, q: f64) -> Result<types::Result, Error> {
        Self::b(n, q)
    }

    /// Store a series of Mathieu characteristic values $a_n(q)$, for
    /// $n$ from `order_min` to `order_max` inclusive, in the array
    /// `result_array`.
    #[doc(alias = "gsl_sf_mathieu_a_array")]
    pub fn a_array(
        &mut self,
        order_min: i32,
        order_max: i32,
        q: f64,
        result_array: &mut [f64],
    ) -> Result<(), Error> {
        let len = order_max - order_min;
        if len < 0 || len as usize > result_array.len() {
            return Err(Error::Invalid);
        }
        let ret = unsafe {
            sys::gsl_sf_mathieu_a_array(
                order_min,
                order_max,
                q,
                self.unwrap_unique(),
                result_array.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }

    /// Store a series of Mathieu characteristic values $b_n(q)$ for
    /// $n$ from `order_min` to `order_max` inclusive, in the array
    /// `result_array`.
    #[doc(alias = "gsl_sf_mathieu_b_array")]
    pub fn b_array(
        &mut self,
        order_min: i32,
        order_max: i32,
        q: f64,
        result_array: &mut [f64],
    ) -> Result<(), Error> {
        let len = order_max - order_min;
        if len < 0 || len as usize > result_array.len() {
            return Err(Error::Invalid);
        }
        let ret = unsafe {
            sys::gsl_sf_mathieu_b_array(
                order_min,
                order_max,
                q,
                self.unwrap_unique(),
                result_array.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }

    /// Return the angular Mathieu function $\ce_n(q,x)$.
    #[doc(alias = "gsl_sf_mathieu_ce_e")]
    pub fn ce(n: i32, q: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_mathieu_ce_e(n, q, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// Return the angular Mathieu function $\se_n(q,x)$.
    #[doc(alias = "gsl_sf_mathieu_se_e")]
    pub fn se(n: i32, q: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_mathieu_se_e(n, q, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// Store a series of the angular Mathieu function $\ce_n(q,x)$ of
    /// order $n$ from `nmin` to `nmax` inclusive, in the array
    /// `result_array`.
    #[doc(alias = "gsl_sf_mathieu_ce_array")]
    pub fn ce_array(
        &mut self,
        nmin: i32,
        nmax: i32,
        q: f64,
        x: f64,
        result_array: &mut [f64],
    ) -> Result<(), Error> {
        let len = nmax - nmin;
        if len < 0 || len as usize > result_array.len() {
            return Err(Error::Invalid);
        }
        let ret = unsafe {
            sys::gsl_sf_mathieu_ce_array(
                nmin,
                nmax,
                q,
                x,
                self.unwrap_unique(),
                result_array.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }

    /// Store a series of the angular Mathieu function and
    /// $\se_n(q,x)$ of order $n$ from `nmin` to `nmax` inclusive, in
    /// the array `result_array`.
    #[doc(alias = "gsl_sf_mathieu_se_array")]
    pub fn se_array(
        &mut self,
        nmin: i32,
        nmax: i32,
        q: f64,
        x: f64,
        result_array: &mut [f64],
    ) -> Result<(), Error> {
        let len = nmax - nmin;
        if len < 0 || len as usize > result_array.len() {
            return Err(Error::Invalid);
        }
        let ret = unsafe {
            sys::gsl_sf_mathieu_se_array(
                nmin,
                nmax,
                q,
                x,
                self.unwrap_unique(),
                result_array.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }

    /// This routine computes the radial `j`-th kind Mathieu function
    /// $\Mc_n^{(j)}(q,x)$ of order `n`.
    ///
    /// The allowed values of `j` are 1 and 2.  The functions for `j`
    /// = 3,4 can be computed as $\Mc_n^{(3)} = \Mc_n^{(1)} + i
    /// \Mc_n^{(2)}$ and $\Mc_n^{(4)} = \Mc_n^{(1)} - i \Mc_n^{(2)}$.
    #[doc(alias = "gsl_sf_mathieu_Mc_e")]
    pub fn Mc(j: i32, n: i32, q: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_mathieu_Mc_e(j, n, q, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the radial `j`-th kind Mathieu function
    /// $\Ms_n^{(j)}(q,x)$ of order `n`.
    ///
    /// The allowed values of `j` are 1 and 2.  The functions for `j`
    /// = 3,4 can be computed as $\Ms_n^{(3)} = \Ms_n^{(1)} + i
    /// \Ms_n^{(2)}$ and $\Ms_n^{(4)} = \Ms_n^{(1)} - i \Ms_n^{(2)}$.
    #[doc(alias = "gsl_sf_mathieu_Ms_e")]
    pub fn Ms(j: i32, n: i32, q: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_mathieu_Ms_e(j, n, q, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// Store a series of the radial Mathieu functions of kind `j`,
    /// with order from `nmin` to `nmax` inclusive, in the array
    /// `result_array`.
    #[doc(alias = "gsl_sf_mathieu_Mc_array")]
    pub fn Mc_array(
        &mut self,
        j: i32,
        nmin: i32,
        nmax: i32,
        q: f64,
        x: f64,
        result_array: &mut [f64],
    ) -> Result<(), Error> {
        let len = nmax - nmin;
        if len < 0 || len as usize > result_array.len() {
            return Err(Error::Invalid);
        }
        let ret = unsafe {
            sys::gsl_sf_mathieu_Mc_array(
                j,
                nmin,
                nmax,
                q,
                x,
                self.unwrap_unique(),
                result_array.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }

    /// Store a series of the radial Mathieu functions of kind `j`,
    /// with order from `nmin` to `nmax` inclusive, in the array
    /// `result_array`.
    #[doc(alias = "gsl_sf_mathieu_Ms_array")]
    pub fn Ms_array(
        &mut self,
        j: i32,
        nmin: i32,
        nmax: i32,
        q: f64,
        x: f64,
        result_array: &mut [f64],
    ) -> Result<(), Error> {
        let len = nmax - nmin;
        if len < 0 || len as usize > result_array.len() {
            return Err(Error::Invalid);
        }
        let ret = unsafe {
            sys::gsl_sf_mathieu_Ms_array(
                j,
                nmin,
                nmax,
                q,
                x,
                self.unwrap_unique(),
                result_array.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }
}
