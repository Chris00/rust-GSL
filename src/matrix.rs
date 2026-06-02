//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Matrices

Matrices are represented by [`MatF32`], [`MatF64`], [`MatI32`], and
[`MatU32`] depending of the type of elements the hold.

## Example

```
use rgsl::MatF64;
let mut m: MatF64 = [[1., 2.], [3., 4.]].into();
assert_eq!(m[(0,0)], 1.);
assert_eq!(m[(0,1)], 2.);
assert_eq!(m[(1,0)], 3.);
assert_eq!(m[(1,1)], 4.);
m[(0,1)] = m[(0,1)] + m[(1,0)];
assert_eq!(m[(0,1)], 5.);
```

 */

use crate::{
    Error,
    ffi::{self, FFI},
    vector::{VecF32, VecF64, VecI32, VecU32},
    view::{View, ViewMut},
};
use pastey::paste;
use std::{
    fmt::{self, Debug, Formatter},
    ops::{Deref, DerefMut, Index, IndexMut},
    ptr,
};

/// Used to describe types accepted as matrices.
///
/// Matrices are stored in row-major order, meaning that each row
/// of elements forms a contiguous block in memory.  This is the
/// standard “C-language ordering” of two-dimensional arrays.
/// Note that FORTRAN stores arrays in column-major order.  The
/// number of rows is `nrows`. The range of valid row indices runs
/// from `0` to `nrows - 1`.  Similarly `ncols` is the number of
/// columns.  The range of valid column indices runs from `0` to
/// `ncols - 1`.  The physical row dimension `tda`, or trailing
/// dimension, specifies the size of a row of the matrix as laid
/// out in memory.
///
/// Matrices are stored in row-major order.
/// For example, in the following matrix `nrows` is 3, `ncols` is
/// 4, and `tda` is 8.  The physical memory layout of the matrix
/// begins in the top left hand-corner and proceeds from left to
/// right along each row in turn.
///
/// ```text
/// 00 01 02 03 XX XX XX XX
/// 10 11 12 13 XX XX XX XX
/// 20 21 22 23 XX XX XX XX
/// ```
///
/// Each unused memory location is represented by “XX”.  The first
/// element of the slice `data` gives the location of the first
/// element of the matrix in memory.
pub trait Matrix<F> {
    /// Number of rows.
    fn nrows(m: &Self) -> usize;
    /// Number of columns.
    fn ncols(m: &Self) -> usize;
    /// Memory slots used for each row.
    fn tda(m: &Self) -> usize;

    /// The memory slice where the matrix is allocated.
    ///
    /// The element (0, 0) of the matrix is at index 0.
    fn as_slice(m: &Self) -> &[F];
}

pub trait MatrixMut<F>: Matrix<F> {
    /// Same as [`Matrix::as_slice`] but for mutale data.
    fn as_mut_slice(m: &mut Self) -> &mut [F];
}

/// Convert `m` to a GSL matrix.
#[allow(dead_code)]
pub(crate) fn matrix_as_gsl<M>(m: &M) -> View<'_, sys::gsl_matrix>
where
    M: Matrix<f64> + ?Sized,
{
    let m = sys::gsl_matrix {
        size1: Matrix::nrows(m),
        size2: Matrix::ncols(m),
        tda: Matrix::tda(m),
        data: Matrix::as_slice(m).as_ptr() as *mut _,
        block: ptr::null_mut(),
        owner: 0,
    };
    View::new(m, false)
}

/// Convert `m` to a GSL matrix.
#[allow(dead_code)]
pub(crate) fn matrix_as_gsl_mut<M>(m: &mut M) -> ViewMut<'_, sys::gsl_matrix>
where
    M: MatrixMut<f64> + ?Sized,
{
    let m = sys::gsl_matrix {
        size1: Matrix::nrows(m),
        size2: Matrix::ncols(m),
        tda: Matrix::tda(m),
        data: MatrixMut::as_mut_slice(m).as_mut_ptr(),
        block: ptr::null_mut(),
        owner: 0,
    };
    ViewMut::new(m, false)
}

/// Used to represent how to convert the GSL matrix format to the
/// matrix (view) type associated with the vector type `Self`.
///
/// It must be implemented on “matrix references” to which other
/// matrix types [`Deref`].
pub trait AsMatrix {
    type MatView<'a>: Deref<Target = Self>;
    type MatViewMut<'a>: DerefMut<Target = Self>;

    /// Return a mutable matrix view (i.e. non-owning) from the slice.
    ///
    /// See [`Matrix`] for more information on the memory layout.
    fn mat_view_from_slice(
        data: &[f64],
        nrows: usize,
        ncols: usize,
        tda: usize,
    ) -> Self::MatView<'_>;

    /// Same as [`Self::mat_view_from_slice`] but for mutable data.
    fn mat_view_from_mut_slice(
        data: &mut [f64],
        nrows: usize,
        ncols: usize,
        tda: usize,
    ) -> Self::MatViewMut<'_>;

    /// Convert the mutable GSL pointer as a mutable matrix view.
    ///
    /// # Safety
    ///
    /// The GSL vector is not owned and neither the C struct nor the
    /// data block must be freed.
    ///
    /// It is important to ensure that the view lifetime is bound to
    /// the lifetime of the vector or matrix that underlies `view`.
    unsafe fn mat_view_from_ptr<'a>(ptr: *const sys::gsl_matrix) -> Self::MatView<'a> {
        debug_assert!(!ptr.is_null());
        let v = unsafe { ptr.as_ref_unchecked() };
        let nrows = v.size1;
        let ncols = v.size2;
        let len = v.tda * nrows;
        let data = unsafe { std::slice::from_raw_parts(v.data, len) };
        Self::mat_view_from_slice(data, nrows, ncols, v.tda)
    }

    /// Same as [`Self::mat_view_from_ptr`] but for mutable data.
    ///
    /// # Safety
    ///
    /// The GSL matrix is not owned and neither the C struct nor the
    /// data block must be freed.
    ///
    /// It is important to ensure that the view lifetime is bound to
    /// the lifetime of the vector or matrix that underlies `view`.
    unsafe fn mat_view_from_mut_ptr<'a>(ptr: *mut sys::gsl_matrix) -> Self::MatViewMut<'a> {
        debug_assert!(!ptr.is_null());
        let v = unsafe { ptr.as_ref_unchecked() };
        let nrows = v.size1;
        let ncols = v.size2;
        let len = v.tda * nrows;
        let data = unsafe { std::slice::from_raw_parts_mut(v.data, len) };
        Self::mat_view_from_mut_slice(data, nrows, ncols, v.tda)
    }
}

impl AsMatrix for MatF64 {
    type MatView<'a> = View<'a, MatF64>;
    type MatViewMut<'a> = ViewMut<'a, MatF64>;

    fn mat_view_from_slice(
        data: &[f64],
        nrows: usize,
        ncols: usize,
        tda: usize,
    ) -> Self::MatView<'_> {
        let m = sys::gsl_matrix {
            size1: nrows,
            size2: ncols,
            tda,
            data: data.as_ptr() as *mut f64, // `View` disallows mutation
            block: ptr::null_mut(),
            owner: 0,
        };
        View::alloc("<MatF64 as AsMatrix>::view_from_slice", m)
    }

    fn mat_view_from_mut_slice(
        data: &mut [f64],
        nrows: usize,
        ncols: usize,
        tda: usize,
    ) -> Self::MatViewMut<'_> {
        let m = sys::gsl_matrix {
            size1: nrows,
            size2: ncols,
            tda,
            data: data.as_mut_ptr(),
            block: ptr::null_mut(),
            owner: 0,
        };
        ViewMut::alloc("<MatF64 as AsMatrix>::view_from_mut_slice", m)
    }

    unsafe fn mat_view_from_ptr<'a>(ptr: *const sys::gsl_matrix) -> Self::MatView<'a> {
        View::from_ptr(ptr, false)
    }

    unsafe fn mat_view_from_mut_ptr<'a>(ptr: *mut sys::gsl_matrix) -> Self::MatViewMut<'a> {
        ViewMut::from_ptr(ptr, false)
    }
}

macro_rules! gsl_matrix {
    ($rust_name:ident, $name:ident, $rust_ty:ident, $vec_name:ident, $vec_c_name:ident) => {
        paste! {
            #[doc = "A matrix with elements of type [`" $rust_ty "`]."]
            pub struct $rust_name {
                mat: *mut sys::$name,
            }

            impl $rust_name {
                #[doc = "Creates a new " $rust_name " with all elements set to zero"]
                #[doc(alias = $name _calloc)]
                pub fn new(nrows: usize, ncols: usize) -> Self {
                    let tmp = unsafe { sys::[<$name _calloc>](nrows, ncols) };
                    if tmp.is_null() {
                        panic!("{}::new cannot allocate memory",
                            stringify!($rust_name));
                    }
                    Self::wrap(tmp)
                }

                /// Convert the GSL view.
                ///
                /// SAFETY: It is important to ensure that the view lifetime
                /// is bound to the lifetime of the vector or matrix that
                /// underlies `view`.
                pub(crate) unsafe fn view_mut<'a>(
                    view: sys::[<$name _view>]
                ) -> ViewMut<'a, Self> {
                    // The view contains a stack allocated `gsl_matrix`.
                    // Since we want to be compatible `MatXXX`, we reallocate
                    // it on the heap.
                    ViewMut::alloc("matrix::view_mut", view.matrix)
                }

                /// Return a matrix view of `base`.  The matrix has
                /// `nrows` rows and `ncols` columns.  The (`i`, `j`)-th
                /// element of the new matrix is given by
                ///
                /// ```text
                /// m'(i,j) = base[i*tda + j]
                /// ```
                ///
                /// where the index `i` runs from `0` to `nrows-1` and
                /// the index `j` runs from `0` to `ncols-1`.
                ///
                /// Panic if `nrows * ncols` is larger than the length
                /// of `base`.
                #[doc(alias = $name _const_view_array)]
                pub fn from_slice(
                    base: &[$rust_ty],
                    nrows: usize,
                    ncols: usize,
                ) -> View<'_, Self> {
                    if nrows * ncols > base.len() {
                        panic!("rgsl::matrix::from_slice: nrows * ncols = {nrows} * {ncols} > slice length = {}",
                            base.len())
                    }
                    // `sys::$name _const_view_array` only creates a
                    // view on the stack.  Since we want a
                    // representation compatible with a GSL vector, we
                    // allocate (on the heap) a vector with a NULL
                    // non-owned block — as GSL view_array does.
                    // These views must drop their value in order for
                    // the gsl_vector struct to be freed.
                    View::alloc("matrix::from_slice",
                        sys::$name {
                            size1: nrows,
                            size2: ncols,
                            tda: ncols,
                            // Mutability change is Ok because `View`
                            // does not allow mutability.
                            data: base.as_ptr() as *mut $rust_ty,
                            block: ptr::null_mut(),
                            owner: 0,
                        })
                }


                #[doc(alias = $name _view_array)]
                pub fn from_mut_slice(
                    base: &mut [$rust_ty],
                    nrows: usize,
                    ncols: usize,
                ) -> ViewMut<'_, Self> {
                    if nrows * ncols > base.len() {
                        panic!("rgsl::matrix::from_mut_slice: nrows * ncols = {nrows} * {ncols} > slice length = {}",
                            base.len())
                    }
                    // `sys::$name _view_array` only creates a view on the
                    // stack.  Since we want a representation compatible with
                    // a GSL vector, we allocate (on the heap) a vector with a
                    // NULL non-owned block — as GSL view_array does.  These
                    // views must drop their value in order for the gsl_vector
                    // struct to be freed.
                    ViewMut::alloc("matrix::from_mut_slice",
                        sys::$name {
                            size1: nrows,
                            size2: ncols,
                            tda: ncols,
                            data: base.as_mut_ptr(),
                            block: ptr::null_mut(),
                            owner: 0,
                        })
                }

                /// Return the (`i`, `j`)-th element of the matrix.
                /// If `i` or `j` lie outside the allowed range of `0`
                /// to `nrows - 1` and `0` to `ncols - 1` respectively
                /// then the error handler is invoked and 0 is
                /// returned.
                #[doc(alias = $name _get)]
                pub fn get(&self, i: usize, j: usize) -> $rust_ty {
                    unsafe { sys::[<$name _get>](self.unwrap_shared(), i, j) }
                }

                /// This function sets the value of the (`i`, `j`)-th
                /// element of the matrix to `value`.  If `i` or `j`
                /// lies outside the allowed range of `0` to `nrows -
                /// 1` and `0` to `ncols - 1` then the error handler
                /// is invoked.
                #[doc(alias = $name _set)]
                pub fn set(&mut self, i: usize, j: usize, value: $rust_ty) -> &$rust_name {
                    unsafe { sys::[<$name _set>](self.unwrap_unique(), i, j, value) };
                    self
                }

                /// Set all the elements of the matrix to the value `x`.
                #[doc(alias = $name _set_all)]
                pub fn set_all(&mut self, x: $rust_ty) -> &$rust_name {
                    unsafe { sys::[<$name _set_all>](self.unwrap_unique(), x) };
                    self
                }

                /// Set all the elements of the matrix to zero.
                #[doc(alias = $name _set_zero)]
                pub fn set_zero(&mut self) -> &$rust_name {
                    unsafe { sys::[<$name _set_zero>](self.unwrap_unique()) };
                    self
                }

                /// Set the elements of the matrix to the
                /// corresponding elements of the identity matrix,
                /// `m(i,j)`$= δᵢⱼ$, i.e. a unit diagonal with all
                /// off-diagonal elements zero.  This applies to both
                /// square and rectangular matrices.
                #[doc(alias = $name _set_identity)]
                pub fn set_identity(&mut self) -> &$rust_name {
                    unsafe { sys::[<$name _set_identity>](self.unwrap_unique()) };
                    self
                }

                /// Copy the elements of the `other` matrix into the
                /// `self` matrix.  The two matrices must have the same size.
                #[doc(alias = $name _memcpy)]
                pub fn copy_from(&mut self, other: &$rust_name) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _memcpy>](
                        self.unwrap_unique(), other.unwrap_shared()) };
                    Error::handle(ret, ())
                }

                /// Copy the elements of the `self` matrix into the
                /// `other` matrix.  The two matrices must have the same size.
                #[doc(alias = $name _memcpy)]
                pub fn copy_to(&self, other: &mut $rust_name) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _memcpy>](
                        other.unwrap_unique(), self.unwrap_shared()) };
                    Error::handle(ret, ())
                }

                /// Exchange the elements of the matrices self and other by
                /// copying.  The two matrices must have the same size.
                #[doc(alias = $name _swap)]
                pub fn swap(&mut self, other: &mut $rust_name) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _swap>](
                        self.unwrap_unique(), other.unwrap_unique()) };
                    Error::handle(ret, ())
                }

                /// Copy the elements of the `i`-th row of the
                /// matrix into the returned vector.
                #[doc(alias = $name _get_row)]
                pub fn get_row(&self, i: usize) -> Result<$vec_name, Error> {
                    let tmp = unsafe { sys::[<$vec_c_name _alloc>](self.ncols()) };

                    if tmp.is_null() {
                        // TODO: So we want this or to panic?
                        Err(Error::NoMemory)
                    } else {
                        let ret = unsafe { sys::[<$name _get_row>](
                            tmp, self.unwrap_shared(), i) };

                        Error::handle(ret, ffi::FFI::wrap(tmp))
                    }
                }

                /// Copy the elements of the `j`-th column of the
                /// matrix into the returned vector.
                #[doc(alias = $name _get_col)]
                pub fn get_col(&self, j: usize) -> Result<$vec_name, Error> {
                    let tmp = unsafe { sys::[<$vec_c_name _alloc>](self.nrows()) };

                    if tmp.is_null() {
                        Err(Error::NoMemory)
                    } else {
                        let ret = unsafe { sys::[<$name _get_col>](
                            tmp, self.unwrap_shared(), j) };

                        Error::handle(ret, ffi::FFI::wrap(tmp))
                    }
                }

                /// Copy the elements of the vector `v` into the
                /// `i`-th row of the matrix.  The length of the
                /// vector must be the same as the length of the row.
                #[doc(alias = $name _set_row)]
                pub fn set_row(
                    &mut self,
                    i: usize,
                    v: &$vec_name,
                ) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _set_row>](
                        self.unwrap_unique(), i, v.unwrap_shared()) };
                    Error::handle(ret, ())
                }

                /// Copy the elements of the vector `v` into the
                /// `j`-th column of the matrix.  The length of the
                /// vector must be the same as the length of the column.
                #[doc(alias = $name _set_col)]
                pub fn set_col(
                    &mut self,
                    j: usize,
                    v: &$vec_name,
                ) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _set_col>](
                        self.unwrap_unique(), j, v.unwrap_shared()) };
                    Error::handle(ret, ())
                }

                /// Exchange the `i1`th and `i2`th rows of the matrix in-place.
                #[doc(alias = $name _swap_rows)]
                pub fn swap_rows(&mut self, i1: usize, i2: usize) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _swap_rows>](
                        self.unwrap_unique(), i1, i2) };
                    Error::handle(ret, ())
                }

                /// Exchange the `j1`th and `j2`th columns of the matrix
                /// in-place.
                #[doc(alias = $name _swap_columns)]
                pub fn swap_columns(
                    &mut self,
                    j1: usize,
                    j2: usize,
                ) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _swap_columns>](
                        self.unwrap_unique(), j1, j2) };
                    Error::handle(ret, ())
                }

                /// Exchange the `i`-th row and `j`-th column of the
                /// matrix in-place.  The matrix must be square for
                /// this operation to be possible.
                #[doc(alias = $name _swap_rowcol)]
                pub fn swap_row_col(&mut self, i: usize, j: usize) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _swap_rowcol>](self.unwrap_unique(), i, j) };
                    Error::handle(ret, ())
                }

                /// Return the transpose of the matrix by copying the
                /// elements into it.  This function works for all
                /// matrices provided that the dimensions of the
                /// matrix dest match the transposed dimensions of the
                /// matrix.
                #[doc(alias = $name _transpose_memcpy)]
                pub fn transpose_memcpy(&self) -> Result<$rust_name, Error> {
                    let dest = unsafe { sys::[<$name _alloc>](
                        self.ncols(), self.nrows()) };

                    if dest.is_null() {
                        Err(Error::NoMemory)
                    } else {
                        let ret = unsafe { sys::[<$name _transpose_memcpy>](
                            dest, self.unwrap_shared()) };

                        Error::handle(ret, $rust_name::wrap(dest))
                    }
                }

                /// Replace the matrix m by its transpose by copying the
                /// elements of the matrix in-place. The matrix must be square for
                /// this operation to be possible.
                #[doc(alias = $name _transpose)]
                pub fn transpose(&mut self) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _transpose>](
                        self.unwrap_unique()) };
                    Error::handle(ret, ())
                }

                /// Add the elements of the `other` matrix to the
                /// elements of the `self` matrix.  The result
                /// `self(i,j) ← self(i,j) + other(i,j)` is stored in
                /// `self` and `other` remains unchanged.  The two
                /// matrices must have the same dimensions.
                #[doc(alias = $name _add)]
                pub fn add(&mut self, other: &$rust_name) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _add>](
                        self.unwrap_unique(),
                        other.unwrap_shared()
                    )};
                    Error::handle(ret, ())
                }

                /// Subtract the elements of the `other` matrix from
                /// the elements of the `self` matrix.  The result
                /// `self(i,j) ← self(i,j) - other(i,j)` is stored in
                /// `self` and `other` remains unchanged.  The two
                /// matrices must have the same dimensions.
                #[doc(alias = $name _sub)]
                pub fn sub(&mut self, other: &$rust_name) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _sub>](
                        self.unwrap_unique(),
                        other.unwrap_shared()
                    )};
                    Error::handle(ret, ())
                }

                /// Multiply the elements of the `self` matrix by the
                /// elements of the `other` matrix.  The result
                /// `self(i,j) ← self(i,j) * other(i,j)` is stored in
                /// `self` and `other` remains unchanged.  The two
                /// matrices must have the same dimensions.
                #[doc(alias = $name _mul_elements)]
                pub fn mul_elements(&mut self, other: &$rust_name) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _mul_elements>](
                        self.unwrap_unique(),
                        other.unwrap_shared(),
                    )};
                    Error::handle(ret, ())
                }

                /// Divide the elements of the `self` matrix by the
                /// elements of the `other` matrix.  The result
                /// `self(i,j) ← self(i,j) / other(i,j)` is stored in
                /// `self` and `other` remains unchanged.  The two
                /// matrices must have the same dimensions.
                #[doc(alias = $name _div_elements)]
                pub fn div_elements(&mut self, other: &$rust_name) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _div_elements>](
                        self.unwrap_unique(),
                        other.unwrap_shared(),
                    )};
                    Error::handle(ret, ())
                }

                /// Multiply the elements of the `self` matrix by the
                /// constant factor `x`.  The result `self(i,j) ← x
                /// self(i,j)` is stored in `self`.
                #[doc(alias = $name _scale)]
                pub fn scale(&mut self, x: $rust_ty) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _scale>](
                        self.unwrap_unique(), x) };
                    Error::handle(ret, ())
                }

                /// Add the constant value `x` to the elements of the
                /// self matrix.  The result `self(i,j) ← self(i,j) + x`
                /// is stored in `self`.
                #[doc(alias = $name _add_constant)]
                pub fn add_constant(&mut self, x: $rust_ty) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _add_constant>](
                        self.unwrap_unique(), x) };
                    Error::handle(ret, ())
                }

                #[doc(alias = $name _add_diagonal)]
                pub fn add_diagonal(&mut self, x: $rust_ty) -> Result<(), Error> {
                    let ret = unsafe { sys::[<$name _add_diagonal>](
                        self.unwrap_unique(), x) };
                    Error::handle(ret, ())
                }

                /// Return the maximum value in the matrix.
                #[doc(alias = $name _max)]
                pub fn max(&self) -> $rust_ty {
                    unsafe { sys::[<$name _max>](self.unwrap_shared()) }
                }

                /// Return the minimum value in the matrix.
                #[doc(alias = $name _min)]
                pub fn min(&self) -> $rust_ty {
                    unsafe { sys::[<$name _min>](self.unwrap_shared()) }
                }

                /// Return the minimum and maximum values in the matrix.
                #[doc(alias = $name _minmax)]
                pub fn minmax(&self) -> ($rust_ty, $rust_ty) {
                    let mut min_out = 0 as _;
                    let mut max_out = 0 as _;
                    unsafe { sys::[<$name _minmax>](
                        self.unwrap_shared(), &mut min_out, &mut max_out) };
                    (min_out, max_out)
                }

                /// Return the indices of the maximum value in the
                /// matrix. When there are several equal maximum
                /// elements then the first element found is returned,
                /// searching in row-major order.
                #[doc(alias = $name _max_index)]
                pub fn max_index(&self) -> (usize, usize) {
                    let mut imax = 0;
                    let mut jmax = 0;
                    unsafe { sys::[<$name _max_index>](
                        self.unwrap_shared(), &mut imax, &mut jmax) };
                    (imax, jmax)
                }

                /// Return the indices of the minimum value in the
                /// matrix.  When there are several equal minimum
                /// elements then the first element found is returned,
                /// searching in row major order.
                #[doc(alias = $name _min_index)]
                pub fn min_index(&self) -> (usize, usize) {
                    let mut imax = 0;
                    let mut jmax = 0;
                    unsafe { sys::[<$name _min_index>](
                        self.unwrap_shared(), &mut imax, &mut jmax) };
                    (imax, jmax)
                }

                /// Return `(imin, jmin, imax, jmax)`, the indices of
                /// the minimum and maximum values in the matrix.
                /// When there are several equal minimum or maximum
                /// elements then the first elements found are
                /// returned, searching in row-major order.
                #[doc(alias = $name _minmax_index)]
                pub fn minmax_index(&self) -> (usize, usize, usize, usize) {
                    let mut imin = 0;
                    let mut jmin = 0;
                    let mut imax = 0;
                    let mut jmax = 0;
                    unsafe { sys::[<$name _minmax_index>](
                        self.unwrap_shared(),
                        &mut imin,
                        &mut jmin,
                        &mut imax,
                        &mut jmax,
                    )};
                    (imin, jmin, imax, jmax)
                }

                /// Return `true` if all the elements of the self
                /// matrix are zero.
                #[doc(alias = $name _isnull)]
                pub fn is_null(&self) -> bool {
                    unsafe { sys::[<$name _isnull>](self.unwrap_shared()) == 1 }
                }

                /// Return `true` if all the elements of the self matrix are
                /// stricly positive.
                #[doc(alias = $name _ispos)]
                pub fn is_pos(&self) -> bool {
                    unsafe { sys::[<$name _ispos>](self.unwrap_shared()) == 1 }
                }

                /// Return `true` if all the elements of the self matrix are
                /// stricly negative.
                #[doc(alias = $name _isneg)]
                pub fn is_neg(&self) -> bool {
                    unsafe { sys::[<$name _isneg>](self.unwrap_shared()) == 1 }
                }

                /// Return `true` if all the elements of the self matrix are
                /// non-negative.
                #[doc(alias = $name _isnonneg)]
                pub fn is_non_neg(&self) -> bool {
                    unsafe { sys::[<$name _isnonneg>](self.unwrap_shared()) == 1 }
                }

                /// Return `true` if all elements of the two matrix are equal.
                #[doc(alias = $name _equal)]
                pub fn equal(&self, other: &$rust_name) -> bool {
                    unsafe {
                        sys::[<$name _equal>](
                            self.unwrap_shared(), other.unwrap_shared()) == 1 }
                }

                #[doc(alias = $name _row)]
                pub fn row(&mut self, i: usize) -> ViewMut<'_, $vec_name> {
                    unsafe {
                        $vec_name::view_mut(sys::[<$name _row>](
                            self.unwrap_unique(), i))
                    }
                }

                #[doc(alias = $name _column)]
                pub fn column(&mut self, j: usize) -> ViewMut<'_, $vec_name> {
                    unsafe {
                        $vec_name::view_mut(sys::[<$name _column>](
                            self.unwrap_unique(), j))
                    }
                }

                #[doc(alias = $name _diagonal)]
                pub fn diagonal(&mut self) -> ViewMut<'_, $vec_name> {
                    unsafe {
                        $vec_name::view_mut(sys::[<$name _diagonal>](
                            self.unwrap_unique()))
                    }
                }

                #[doc(alias = $name _subdiagonal)]
                pub fn subdiagonal(&mut self, k: usize) -> ViewMut<'_, $vec_name> {
                    unsafe {
                        $vec_name::view_mut(sys::[<$name _subdiagonal>](
                            self.unwrap_unique(), k))
                    }
                }

                #[doc(alias = $name _superdiagonal)]
                pub fn superdiagonal(&mut self, k: usize) -> ViewMut<'_, $vec_name> {
                    unsafe {
                        $vec_name::view_mut(sys::[<$name _superdiagonal>](
                            self.unwrap_unique(), k))
                    }
            }

                #[doc(alias = $name _subrow)]
                pub fn subrow(
                    &mut self,
                    i: usize,
                    offset: usize,
                    n: usize,
                ) -> ViewMut<'_, $vec_name> {
                    unsafe {
                        $vec_name::view_mut(sys::[<$name _subrow>](
                            self.unwrap_unique(), i, offset, n))
                    }
                }

                #[doc(alias = $name _subcolumn)]
                pub fn subcolumn(
                    &mut self,
                    i: usize,
                    offset: usize,
                    n: usize,
                ) -> ViewMut<'_, $vec_name> {
                    unsafe {
                        $vec_name::view_mut(sys::[<$name _subcolumn>](
                            self.unwrap_unique(), i, offset, n))
                    }
                }

                /// Return a matrix view of a submatrix.  The
                /// upper-left element of the submatrix is the element
                /// (`k1`, `k2`) of the original matrix.  The
                /// submatrix has `nrows` rows and `ncols` columns.
                ///
                /// Panic if the combined parameters (`k1`, `k2`,
                /// `nrows`, `ncols`) overrun the ends of the original
                /// matrix.
                ///
                // TODO: The function gsl_matrix_const_submatrix is equivalent
                // to gsl_matrix_submatrix but can be used for matrices which
                // are declared const.
                #[doc(alias = $name _submatrix)]
                pub fn submatrix(
                    &mut self,
                    k1: usize,
                    k2: usize,
                    nrows: usize,
                    ncols: usize,
                ) -> ViewMut<'_, Self> {
                    unsafe {
                        let view = sys::[<$name _submatrix>](self.mat, k1, k2, nrows, ncols);
                        if view.matrix.data.is_null() {
                            panic!("rgls::matrix::submatrix: (k1, k2, nrows, ncols) = ({k1}, {k2}, {nrows}, {ncols}) exceeds the matrix size ({}, {})",
                                self.nrows(), self.ncols())
                        }
                        $rust_name::view_mut(view)
                    }
                }

                #[inline]
                pub fn nrows(&self) -> usize {
                    debug_assert!(!self.unwrap_shared().is_null());
                    unsafe { (*self.unwrap_shared()).size1 }
                }

                #[inline]
                pub fn ncols(&self) -> usize {
                    debug_assert!(!self.unwrap_shared().is_null());
                    unsafe { (*self.unwrap_shared()).size2 }
                }

                #[inline]
                pub fn tda(&self) -> usize {
                    debug_assert!(!self.unwrap_shared().is_null());
                    unsafe { *self.unwrap_shared() }.tda
                }

                pub fn as_slice(&self) -> &[$rust_ty] {
                    let m = unsafe { *self.unwrap_shared() };
                    let len = m.size1 * m.tda;
                    unsafe { std::slice::from_raw_parts(m.data, len) }
                }

                pub fn as_mut_slice(&mut self) -> &mut [$rust_ty] {
                    let m = unsafe { *self.unwrap_unique() };
                    let len = m.size1 * m.tda;
                    unsafe { std::slice::from_raw_parts_mut(m.data, len) }
                }


                // TODO: impl Clone
                pub fn clone(&self) -> Result<Self, Error> {
                    debug_assert!(!self.unwrap_shared().is_null());
                    let mut m = Self::new(self.nrows(), self.ncols());
                    m.copy_from(self)?;
                    Ok(m)
                }

                // FIXME: why public?
                #[doc(hidden)]
                pub fn is_ptr_null(&self) -> bool {
                    self.unwrap_shared().is_null()
                }
            }

            impl Drop for $rust_name {
                #[doc(alias = $name _free)]
                fn drop(&mut self) {
                    unsafe { sys::[<$name _free>](self.mat) };
                }
            }

            unsafe impl FFI for $rust_name {
                type Sys = sys::$name;

                fn wrap(mat: *mut sys::$name) -> Self {
                    Self { mat }
                }

                fn unwrap_shared(&self) -> *const sys::$name {
                    self.mat as *const _
                }

                fn unwrap_unique(&mut self) -> *mut sys::$name {
                    self.mat
                }
            }

            impl Debug for $rust_name {
                #[allow(unused_must_use)]
                fn fmt(&self, f: &mut Formatter) -> fmt::Result {
                    let ptr = self.unwrap_shared();
                    if ptr.is_null() {
                        write!(f, "<null>")
                    } else {
                        let size1 = self.nrows();
                        let size2 = self.ncols();
                        for y in 0..size1 {
                            write!(f, "[");
                            for x in 0..size2 {
                                if x < size2 - 1 {
                                    write!(f, "{}, ", self.get(y, x));
                                } else {
                                    write!(f, "{}", self.get(y, x));
                                }
                            }
                            if y < size1 - 1 {
                                write!(f, "]\n");
                            }
                        }
                        write!(f, "]")
                    }
                }
            }

            impl Matrix<$rust_ty> for $rust_name {
                fn nrows(m: &Self) -> usize {
                    m.nrows()
                }

                fn ncols(m: &Self) -> usize {
                    m.ncols()
                }

                fn tda(m: &Self) -> usize {
                    m.tda()
                }

                fn as_slice(m: &Self) -> &[$rust_ty] {
                    m.as_slice()
                }

            }

            impl Index<(usize, usize)> for $rust_name {
                type Output = $rust_ty;

                fn index(&self, (i,j): (usize, usize)) -> &$rust_ty {
                    if i >= self.nrows() {
                        panic!("rgsl::matrix::{}: row {} >= {}",
                            stringify!($rust_name), i, self.nrows())
                    }
                    if j >= self.ncols() {
                        panic!("rgsl::matrix::{}: col {} >= {}",
                            stringify!($rust_name), j, self.ncols())
                    }
                    // Implemented "by hand" to return a reference.
                    unsafe{
                        let m: sys::$name = *self.unwrap_shared();
                        &*m.data.add(i * m.tda + j)
                    }
                }
            }

            impl IndexMut<(usize, usize)> for $rust_name {
                fn index_mut(&mut self, (i,j): (usize, usize)) -> &mut $rust_ty {
                    if i >= self.nrows() {
                        panic!("rgsl::matrix::{}: row {} >= {}",
                            stringify!($rust_name), i, self.nrows())
                    }
                    if j >= self.ncols() {
                        panic!("rgsl::matrix::{}: col {} >= {}",
                            stringify!($rust_name), j, self.ncols())
                    }
                    // Implemented "by hand" to return a reference.
                    unsafe{
                        let m: sys::$name = *self.unwrap_unique();
                        &mut *m.data.add(i * m.tda + j)
                    }
                }
            }

            impl<const N: usize, const M: usize> From<[[$rust_ty; M]; N]>
            for $rust_name {
                fn from(value: [[$rust_ty; M]; N]) -> Self {
                    let mut m = Self::new(N, M);
                    for i in 0..m.nrows() {
                        for j in 0..m.ncols() {
                            m[(i, j)] = value[i][j];
                        }
                    }
                    m
                }
            }
        } // end of paste! block
    }; // end of the gsl_matrix macro
}

gsl_matrix!(MatF32, gsl_matrix_float, f32, VecF32, gsl_vector_float);
gsl_matrix!(MatF64, gsl_matrix, f64, VecF64, gsl_vector);
gsl_matrix!(MatI32, gsl_matrix_int, i32, VecI32, gsl_vector_int);
gsl_matrix!(MatU32, gsl_matrix_uint, u32, VecU32, gsl_vector_uint);

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn out_of_bounds() {
        let mut s = [1.];
        let m = MatF64::from_mut_slice(&mut s, 1, 1);
        assert_eq!(m.get(0, 0), 1.);
        assert_eq!(m.get(1, 0), 0.);
    }

    #[test]
    fn indices() {
        let mut s = [1., 2., 3., 4.];
        let mut m = MatF64::from_mut_slice(&mut s, 2, 2);
        assert_eq!(m[(0, 0)], 1.);
        assert_eq!(m[(0, 1)], 2.);
        assert_eq!(m[(1, 0)], 3.);
        assert_eq!(m[(1, 1)], 4.);
        m[(0, 1)] = 5.;
        assert_eq!(m[(0, 1)], 5.);
        drop(m);
        assert_eq!(s[1], 5.);
    }
}
