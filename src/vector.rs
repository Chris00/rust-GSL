//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Vectors

Vectors (resp. mutable vectors) are any type that implements the
[`Vector`] (resp. [`VectorMut`]) trait which defines the interface to a
slice of a memory block.  A vector slice is a set of equally-spaced
elements of an area of memory.
!*/

use crate::{Error, ffi::FFI, view::ViewMut};
use std::fmt::{self, Debug, Formatter};

use pastey::paste;

#[cfg(feature = "complex")]
extern crate num_complex;
#[cfg(feature = "complex")]
use self::num_complex::Complex;

#[allow(clippy::len_without_is_empty)]
/// Trait implemented by types that are considered vectors by this crate.
/// Elements of the vector are of type `F` (e.g., `f32` or `f64`).
///
/// Bring this trait into scope in order to add methods to specify
/// strides to the types implementing `Vector`.
///
/// # Safety
/// One must make sore that `(len - 1) * stride` does not exceed the
/// length of the underlying slice.
pub unsafe trait Vector<F> {
    /// Return the number of elements in the vector.
    ///
    /// This is an associated function rather than a method, in order
    /// to avoid conflicts with methods with the same name.
    fn len(x: &Self) -> usize;

    /// The distance in the slice between two consecutive elements of
    /// the vector in [`Vector::as_slice`] and [`VectorMut::as_mut_slice`].
    fn stride(x: &Self) -> usize;

    /// Return a reference to the underlying slice.  The `i`th element
    /// of the vector, `0 <= i < len(x)`, is the `i * stride` element
    /// in the slice.
    ///
    /// # Safety
    ///
    /// The implementation must ensure that the slice is large enough
    /// to hold all elements of the vector.
    ///
    /// # Remark
    ///
    /// This is an associated function rather than a method in order
    /// to avoid conflicts with methods with the same name.
    fn as_slice(x: &Self) -> &[F];
}

/// Trait implemented by types that are considered *mutable* vectors
/// by this crate.  Elements of the vector are of type `F`
/// (e.g. `f32` or `f64`).
///
/// Bring this trait into scope in order to add methods to specify
/// strides to the types implementing `Vector`.
///
/// # Safety
/// One must make sure that `(len - 1) * stride` does not exceed the
/// length of the underlying slice.
pub unsafe trait VectorMut<F>: Vector<F> {
    /// Same as [`Vector::as_slice`] but mutable.
    fn as_mut_slice(x: &mut Self) -> &mut [F];
}

// Implement the `Vector` trait on standard vectors.

macro_rules! impl_AsRef {
    ($ty: ty) => {
        unsafe impl<T> Vector<$ty> for T
        where
            T: AsRef<[$ty]> + ?Sized,
        {
            #[inline]
            fn len(x: &Self) -> usize {
                x.as_ref().len()
            }
            #[inline]
            fn stride(_: &Self) -> usize {
                1
            }
            #[inline]
            fn as_slice(x: &Self) -> &[$ty] {
                x.as_ref()
            }
        }

        unsafe impl<T> VectorMut<$ty> for T
        where
            T: Vector<$ty> + AsMut<[$ty]> + ?Sized,
        {
            #[inline]
            fn as_mut_slice(x: &mut Self) -> &mut [$ty] {
                x.as_mut()
            }
        }
    };
}

impl_AsRef!(f32);
impl_AsRef!(f64);
#[cfg(feature = "complex")]
impl_AsRef!(Complex<f32>);
#[cfg(feature = "complex")]
impl_AsRef!(Complex<f64>);

/// Return the length of `x` as a `i32` value (to use in CBLAS calls).
#[inline]
pub(crate) fn len<F, T: Vector<F> + ?Sized>(x: &T) -> i32 {
    T::len(x).try_into().expect("Length must fit in `i32`")
}

#[inline]
pub(crate) fn as_ptr<F, T: Vector<F> + ?Sized>(x: &T) -> *const F {
    T::as_slice(x).as_ptr()
}

#[inline]
pub(crate) fn as_mut_ptr<F, T: VectorMut<F> + ?Sized>(x: &mut T) -> *mut F {
    T::as_mut_slice(x).as_mut_ptr()
}

/// Return the stride of `x` as a `i32` value (to use in CBLAS calls).
#[inline]
pub(crate) fn stride<F, T: Vector<F> + ?Sized>(x: &T) -> i32 {
    T::stride(x).try_into().expect("Stride must fit in `i32`")
}

#[inline]
pub(crate) fn check_equal_len<T1, T2, F>(x: &T1, y: &T2) -> Result<(), Error>
where
    T1: Vector<F> + ?Sized,
    T2: Vector<F> + ?Sized,
{
    if T1::len(x) != T2::len(y) {
        return Err(Error::Invalid);
    }
    Ok(())
}

macro_rules! gsl_vec {
    ($rust_name:ident, $name:ident, $rust_ty:ident) => {
paste! {

    #[doc = "Vector with elements of type [`" $rust_ty "`]."]
    pub struct $rust_name {
        vec: *mut sys::$name,
    }

    impl Drop for $rust_name {
        #[doc(alias = $name _free)]
        fn drop(&mut self) {
            unsafe { sys::[<$name _free>](self.unwrap_unique()) };
        }
    }

    impl Debug for $rust_name {
        fn fmt(&self, f: &mut Formatter) -> fmt::Result {
            debug_assert!(! self.unwrap_shared().is_null());
            write!(f, "{:?}", self.as_slice())
        }
    }

    unsafe impl FFI for $rust_name {
        type Sys = sys::$name;

        fn wrap(vec: *mut sys::$name) -> Self {
            Self { vec }
        }

        fn unwrap_shared(&self) -> *const sys::$name {
            self.vec as *const _
        }

        fn unwrap_unique(&mut self) -> *mut sys::$name {
            self.vec
        }
    }

    impl $rust_name {
        #[doc = "Create a new `" $rust_name "` with all elements set to zero."]
        ///
        /// Panic if no memory can be allocated.
        #[doc(alias = $name _calloc)]
        pub fn new(size: usize) -> Self {
            let tmp = unsafe { sys::[<$name _calloc>](size) };
            if tmp.is_null() {
                panic!("{}::new cannot allocate memory",
                    stringify!($rust_name));
            }
            Self::wrap(tmp)
        }

        #[doc = "Create a new [`" $rust_name
            "`] with elements initialized from `slice`."]
        ///
        /// Panic if no memory can be allocated.
        #[doc(alias = $name _alloc)]
        pub fn from_slice(slice: &[$rust_ty]) -> Self {
            let tmp = unsafe { sys::[<$name _alloc>](slice.len() as _) };
            if tmp.is_null() {
                panic!("{}::from_slice cannot allocate memory",
                    stringify!($rust_name));
            }
            let mut v = Self::wrap(tmp);

            for (pos, tmp) in slice.iter().enumerate() {
                v.set(pos as _, *tmp);
            }
            v
        }

        /// Convert the GSL view.
        ///
        /// SAFETY: It is important to ensure that the view lifetime
        /// is bound to the lifetime of the vector or matrix that
        /// underlies `view`.
        pub(crate) unsafe fn view_mut<'a>(
            view: sys::[<$name _view>]
        ) -> ViewMut<'a, Self> {
            // The view contains a stack allocated `gsl_vector`.
            // Since we want to Deref to `VecXXX`, we reallocate
            // it on the heap.
            ViewMut::alloc("vector::view_mut", view.vector)
        }

        pub fn len(&self) -> usize {
            let ptr = self.unwrap_shared();
            debug_assert!(! ptr.is_null());
            unsafe { (*ptr).size }
        }

        pub fn is_empty(&self) -> bool {
            self.len() == 0
        }

        pub fn as_slice(&self) -> &[$rust_ty] {
            let ptr = unsafe { (*self.unwrap_shared()).data };
            debug_assert!(! ptr.is_null());
            unsafe { std::slice::from_raw_parts(ptr, self.len()) }
        }

        pub fn as_slice_mut(&mut self) -> &mut [$rust_ty] {
            let ptr = unsafe { (*self.unwrap_shared()).data };
            debug_assert!(! ptr.is_null());
            unsafe { std::slice::from_raw_parts_mut(ptr, self.len()) }
        }

        /// Return the `i`-th element of the vector `self`.  If `i` lies
        /// outside the allowed range of `0` to `n-1` then 0 is returned.
        #[doc(alias = $name _get)]
        pub fn get(&self, i: usize) -> $rust_ty {
            unsafe { sys::[<$name _get>](self.unwrap_shared(), i) }
        }

        /// Set the value of the `i`-th element of a vector `self` to
        /// `x`.  If `i` lies outside the allowed range of 0 to n-1,
        /// the function panics.
        #[doc(alias = $name _set)]
        pub fn set(&mut self, i: usize, x: $rust_ty) -> &mut $rust_name {
            // TODO: panic
            unsafe { sys::[<$name _set>](self.unwrap_unique(), i, x) };
            self
        }

        /// Set all the elements of the vector `self` to the value `x`.
        #[doc(alias = $name _set_all)]
        pub fn set_all(&mut self, x: $rust_ty) -> &mut $rust_name {
            unsafe { sys::[<$name _set_all>](self.unwrap_unique(), x) };
            self
        }

        /// Set all the elements of the vector `self` to zero.
        #[doc(alias = $name _set_zero)]
        pub fn set_zero(&mut self) -> &mut $rust_name {
            unsafe { sys::[<$name _set_zero>](self.unwrap_unique()) };
            self
        }

        /// Make a basis vector by setting all the elements of the
        /// vector `self` to zero except for the `i`-th element which
        /// is set to one.
        #[doc(alias = $name _set_basis)]
        pub fn set_basis(&mut self, i: usize) -> &mut $rust_name {
            unsafe { sys::[<$name _set_basis>](self.unwrap_unique(), i) };
            self
        }

        /// Copy the elements of the other vector into the `self`
        /// vector.  The two vectors must have the same length.
        #[doc(alias = $name _memcpy)]
        pub fn copy_from(&mut self, other: &$rust_name) -> Result<(), Error> {
            let ret = unsafe { sys::[<$name _memcpy>](
                self.unwrap_unique(),
                other.unwrap_shared()) };
            Error::handle(ret, ())
        }

        /// Copy the elements of the `self` vector into the other
        /// vector.  The two vectors must have the same length.
        #[doc(alias = $name _memcpy)]
        pub fn copy_to(&self, other: &mut $rust_name) -> Result<(), Error> {
            let ret = unsafe { sys::[<$name _memcpy>](
                other.unwrap_unique(),
                self.unwrap_shared()
            )};
            Error::handle(ret, ())
        }

        /// Exchange the elements of the vectors by copying.  The two
        /// vectors must have the same length.
        #[doc(alias = $name _swap)]
        pub fn swap(&mut self, other: &mut $rust_name) -> Result<(), Error> {
            let ret = unsafe { sys::[<$name _swap>](
                other.unwrap_unique(),
                self.unwrap_unique()
            )};
            Error::handle(ret, ())
        }

        /// Exchange the `i`-th and `j`-th elements of the vector
        /// `self` in-place.
        #[doc(alias = $name _swap_elements)]
        pub fn swap_elements(&mut self, i: usize, j: usize) -> Result<(), Error> {
            let ret = unsafe { sys::[<$name _swap_elements>](
                self.unwrap_unique(), i, j) };
            Error::handle(ret, ())
        }

        /// Reverse the order of the elements of the vector.
        #[doc(alias = $name _reverse)]
        pub fn reverse(&mut self) -> Result<(), Error> {
            let ret = unsafe { sys::[<$name _reverse>](self.unwrap_unique()) };
            Error::handle(ret, ())
        }

        /// Add the elements of the other vector to the elements of
        /// the `self` vector.  The result $a_i ← a_i + b_i$ is stored
        /// in `self` and other remains unchanged.  The two vectors
        /// must have the same length.
        #[doc(alias = $name _add)]
        pub fn add(&mut self, other: &$rust_name) -> Result<(), Error> {
            let ret = unsafe { sys::[<$name _add>](
                self.unwrap_unique(),
                other.unwrap_shared()
            )};
            Error::handle(ret, ())
        }

        /// Subtract the elements of the `self` vector from the
        /// elements of the other vector. The result $a_i ← a_i - b_i$
        /// is stored in self and other remains unchanged.  The two
        /// vectors must have the same length.
        #[doc(alias = $name _sub)]
        pub fn sub(&mut self, other: &$rust_name) -> Result<(), Error> {
            let ret = unsafe { sys::[<$name _sub>](
                self.unwrap_unique(),
                other.unwrap_shared()
            )};
            Error::handle(ret, ())
        }

        /// Multiply the elements of the `self` vector a by the
        /// elements of the other vector. The result $a_i ← a_i * b_i$
        /// is stored in `self` and other remains unchanged.  The two
        /// vectors must have the same length.
        #[doc(alias = $name _mul)]
        pub fn mul(&mut self, other: &$rust_name) -> Result<(), Error> {
            let ret = unsafe { sys::[<$name _mul>](
                self.unwrap_unique(),
                other.unwrap_shared()
            ) };
            Error::handle(ret, ())
        }

        /// Divide the elements of the `self` vector by the elements
        /// of the other vector.  The result $a_i ← a_i / b_i$ is
        /// stored in `self` and other remains unchanged.  The two
        /// vectors must have the same length.
        #[doc(alias = $name _div)]
        pub fn div(&mut self, other: &$rust_name) -> Result<(), Error> {
            let ret = unsafe { sys::[<$name _div>](
                self.unwrap_unique(),
                other.unwrap_shared()
            ) };
            Error::handle(ret, ())
        }

        /// Multiply the elements of the `self` vector by the constant
        /// factor `x`.  The result $a_i ← a_i * x$ is stored in `self`.
        #[doc(alias = $name _scale)]
        pub fn scale(&mut self, x: $rust_ty) -> Result<(), Error> {
            let ret = unsafe { sys::[<$name _scale>](self.unwrap_unique(), x) };
            Error::handle(ret, ())
        }

        /// Add the constant value `x` to the elements of the `self`
        /// vector.  The result $a_i ← a_i + x$ is stored in `self`.
        #[doc(alias = $name _add_constant)]
        pub fn add_constant(&mut self, x: $rust_ty) -> Result<(), Error> {
            let ret = unsafe { sys::[<$name _add_constant>](self.unwrap_unique(), x) };
            Error::handle(ret, ())
        }

        /// Return the maximum value in the vector.
        #[doc(alias = $name _max)]
        pub fn max(&self) -> $rust_ty {
            unsafe { sys::[<$name _max>](self.unwrap_shared()) }
        }

        /// Return the minimum value in the self vector.
        #[doc(alias = $name _min)]
        pub fn min(&self) -> $rust_ty {
            unsafe { sys::[<$name _min>](self.unwrap_shared()) }
        }

        /// Return the minimum and maximum values in the self vector.
        #[doc(alias = $name _minmax)]
        pub fn minmax(&self) -> ($rust_ty, $rust_ty) {
            let mut min_out = 0 as _;
            let mut max_out = 0 as _;

            unsafe {
                sys::[<$name _minmax>](
                    self.unwrap_shared(), &mut min_out, &mut max_out);
            }
            (min_out, max_out)
        }

        /// Return the index of the maximum value in the vector.  When
        /// there are several equal maximum elements then the lowest index
        /// is returned.
        #[doc(alias = $name _max_index)]
        pub fn max_index(&self) -> usize {
            unsafe { sys::[<$name _max_index>](self.unwrap_shared()) }
        }

        /// Return the index of the minimum value in the vector.  When
        /// there are several equal minimum elements then the lowest
        /// index is returned.
        #[doc(alias = $name _min_index)]
        pub fn min_index(&self) -> usize {
            unsafe { sys::[<$name _min_index>](self.unwrap_shared()) }
        }

        /// Return the indices of the minimum and maximum values in
        /// the vector.  When there are several equal minimum or
        /// maximum elements then the lowest indices are returned.
        #[doc(alias = $name _minmax_index)]
        pub fn minmax_index(&self) -> (usize, usize) {
            let mut imin = 0;
            let mut imax = 0;

            unsafe { sys::[<$name _minmax_index>](
                self.unwrap_shared(), &mut imin, &mut imax) };
            (imin, imax)
        }

        /// Return `true` if all the elements of the self vector are equal to 0.
        #[doc(alias = $name _isnull)]
        pub fn is_null(&self) -> bool {
            unsafe { sys::[<$name _isnull>](self.unwrap_shared()) == 1 }
        }

        /// Return `true` if all the elements of the self vector are
        /// stricly positive.
        #[doc(alias = $name _ispos)]
        pub fn is_pos(&self) -> bool {
            unsafe { sys::[<$name _ispos>](self.unwrap_shared()) == 1 }
        }

        /// Return `true` if all the elements of the self vector are
        /// stricly negative.
        #[doc(alias = $name _isneg)]
        pub fn is_neg(&self) -> bool {
            unsafe { sys::[<$name _isneg>](self.unwrap_shared()) == 1 }
        }

        /// Return true if all the elements of the self vector are
        /// non-negative.
        #[doc(alias = $name _isnonneg)]
        pub fn is_non_neg(&self) -> bool {
            unsafe { sys::[<$name _isnonneg>](self.unwrap_shared()) == 1 }
        }

        #[doc(alias = $name _equal)]
        pub fn equal(&self, other: &$rust_name) -> bool {
            unsafe { sys::[<$name _equal>](self.unwrap_shared(), other.unwrap_shared()) == 1 }
        }

        // TODO: impl Clone ?
        pub fn clone(&self) -> $rust_name {
            if self.unwrap_shared().is_null() {
                panic!("{}::clone cannot allocate memory",
                    stringify!($rust_name));
            }
            let mut v = $rust_name::new(self.len());
            if v.copy_from(self).is_err() {
                panic!("{}::clone problem copying the vector",
                    stringify!($rust_name));
            }
            v
        }

        // TODO: impl Slice and SliceMut
        #[doc(alias = $name _subvector)]
        pub fn subvector(&mut self, offset: usize, n: usize) -> ViewMut<'_, Self> {
            unsafe {
                let view = sys::[<$name _subvector>](self.vec, offset, n);
                if view.vector.data.is_null() {
                    panic!("rgsl::vector::subvector: (offset, n) = ({offset}, {n}) exceeds the vector length {}", self.len())
                }
                Self::view_mut(view)
            }
        }
    }

    unsafe impl Vector<$rust_ty> for $rust_name {
        #[inline]
        fn len(x: &Self) -> usize {
            $rust_name::len(x)
        }
        #[inline]
        fn stride(x: &Self) -> usize {
            let ptr = x.unwrap_shared();
            if ptr.is_null() {
                1
            } else {
                unsafe { (*ptr).stride }
            }
        }
        #[inline]
        fn as_slice(x: &Self) -> &[$rust_ty] {
            $rust_name::as_slice(x)
        }
    }
    unsafe impl VectorMut<$rust_ty> for $rust_name {
        #[inline]
        fn as_mut_slice(x: &mut Self) -> &mut [$rust_ty] {
            $rust_name::as_slice_mut(x)
        }
    }

} // end of paste! block
}; // end of gsl_vec macro
}

gsl_vec!(VecF32, gsl_vector_float, f32);
gsl_vec!(VecF64, gsl_vector, f64);
gsl_vec!(VecI32, gsl_vector_int, i32);
gsl_vec!(VecU32, gsl_vector_uint, u32);

// https://doc.rust-lang.org/std/primitive.f128.html
//gsl_vec!(VecU128, gsl_vector_long_double, f128);

// Helper trait to convert Complex slices.
pub(crate) trait ComplexSlice<T> {
    // fn as_ptr_fXX(&self) -> *const T;
    fn as_mut_ptr_fXX(&mut self) -> *mut T;
}

macro_rules! impl_ComplexSlice {
    ($ty: ty) => {
        impl ComplexSlice<$ty> for [Complex<$ty>] {
            // fn as_ptr_fXX(&self) -> *const $ty {
            //     // Complex<f64> layout is two consecutive f64 values.
            //     self.as_ptr() as *const $ty
            // }

            fn as_mut_ptr_fXX(&mut self) -> *mut $ty {
                self.as_mut_ptr() as *mut $ty
            }
        }
    };
}

#[cfg(feature = "complex")]
impl_ComplexSlice!(f64);
#[cfg(feature = "complex")]
impl_ComplexSlice!(f32);

#[cfg(test)]
mod test {
    use super::VecF64;

    #[test]
    fn vec_get_out_of_range() {
        let v = VecF64::from_slice(&[1., 2.]);
        assert_eq!(v.get(0), 1.);
        assert_eq!(v.get(2), 0.);
    }

    #[test]
    fn vec_set_out_of_range() {
        let mut v = VecF64::from_slice(&[1., 2.]);
        v.set(2, 0.);
    }
}
