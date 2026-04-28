//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::ffi::FFI;
use std::{
    fmt::{self, Debug, Formatter},
    marker::PhantomData,
    mem::{self, ManuallyDrop},
    ops::{Deref, DerefMut},
};

/// Used to express that the type offers views from `gsl_vector`.
///
/// # Safety
///
/// Make sure the `View` and `ViewMut`, when dropped, do not touch the
/// pointer `Sys`.
pub unsafe trait AsView {
    type Sys;
    type View<'a>
    where
        Self: 'a;
    type ViewMut<'a>
    where
        Self: 'a;

    /// Return a view of `ptr`.
    //fn view(len: usize, stride: usize, data: &[T]) -> Self::View<'_>;
    fn as_view<'a>(ptr: *const Self::Sys) -> Self::View<'a>
    where
        Self: 'a;

    /// Return a mutable view of `ptr`.
    fn as_view_mut<'a>(ptr: *mut Self::Sys) -> Self::ViewMut<'a>
    where
        Self: 'a;
}

// View and ViewMut will be used to
// - wrap raw pointers to gls_vector passed to callbacks;
// - wrap pointers to internal gsl_vectors — e.g. matrix rows —
//   (similar to the first case);
// - hold slices of vector and arrays.

unsafe impl<T: FFI> AsView for T {
    type Sys = <T as FFI>::Sys;
    type View<'a>
        = View<'a, T>
    where
        T: 'a;
    type ViewMut<'a>
        = ViewMut<'a, T>
    where
        T: 'a;

    fn as_view<'a>(ptr: *const Self::Sys) -> Self::View<'a>
    where
        T: 'a,
    {
        View::new(ptr, false)
    }

    fn as_view_mut<'a>(ptr: *mut Self::Sys) -> Self::ViewMut<'a>
    where
        T: 'a,
    {
        ViewMut::new(ptr, false)
    }
}

/// An immutable view to `T`.
pub struct View<'a, T: FFI> {
    // `T` will wrap a pointer.
    inner: ManuallyDrop<T>,
    // - In case of callbacks, the given `gsl_vector...` pointer must
    //   not be freed (so the Rust destructor must not be called).
    // - In case of a slice, the heap allocated `gsl_vector` struct
    //   need to be freed (but it will declare not owning memory it
    //   points to so the vector data will still be valid).
    must_drop: bool,
    phantom: PhantomData<&'a T>,
}

impl<T: FFI> Drop for View<'_, T> {
    fn drop(&mut self) {
        if self.must_drop {
            unsafe { ManuallyDrop::drop(&mut self.inner) };
        }
    }
}

impl<T: FFI> Deref for View<'_, T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T: FFI + Debug> Debug for View<'_, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        <T as Debug>::fmt(self, f)
    }
}

impl<'a, T: FFI + 'a> View<'a, T> {
    pub(crate) fn new(ptr: *const T::Sys, must_drop: bool) -> Self {
        // SAFETY: The conversion to a mutable pointer is fine because
        // `View` will not allow mutable access.
        let t = FFI::wrap(ptr as *mut T::Sys);
        Self {
            inner: ManuallyDrop::new(t),
            must_drop,
            phantom: PhantomData,
        }
    }
}

/// A mutable view to `T`.
pub struct ViewMut<'a, T: FFI> {
    inner: ManuallyDrop<T>,
    must_drop: bool,
    phantom: PhantomData<&'a ()>,
}

impl<T: FFI> Drop for ViewMut<'_, T> {
    fn drop(&mut self) {
        if self.must_drop {
            unsafe { ManuallyDrop::drop(&mut self.inner) };
        }
    }
}

impl<T: FFI> Deref for ViewMut<'_, T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T: FFI> DerefMut for ViewMut<'_, T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}

impl<T: FFI + Debug> Debug for ViewMut<'_, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        <T as Debug>::fmt(self, f)
    }
}

impl<'a, T: FFI + 'a> ViewMut<'a, T> {
    pub(crate) fn new(ptr: *mut T::Sys, must_drop: bool) -> Self {
        let t = FFI::wrap(ptr);
        Self {
            inner: ManuallyDrop::new(t),
            must_drop,
            phantom: PhantomData,
        }
    }

    /// `malloc` a C-struct and fill it with `c`.  Beware that the
    /// value of type `T` will be dropped when the returned value is
    /// so be careful about the resources it owns.
    pub(crate) fn alloc(fn_name: &str, c: T::Sys) -> Self {
        let size = mem::size_of::<T::Sys>();
        let ptr = unsafe { libc::malloc(size) as *mut T::Sys };
        if ptr.is_null() {
            panic!("rgsl::{}: cannot allocate", fn_name);
        }
        unsafe { ptr.write(c) }
        // Need to drop the value for the C-struct to be freed.
        Self::new(ptr, true)
    }
}
