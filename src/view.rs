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

/// An immutable view to `T`.
pub struct View<'a, T> {
    // `T` will wrap a pointer or be a raw GSL struct.
    inner: ManuallyDrop<T>,
    // - In case of callbacks, the given `gsl_vector...` pointer must
    //   not be freed (so the Rust destructor must not be called).
    // - In case of a slice, the heap allocated `gsl_vector` struct
    //   need to be freed (but it will declare not owning memory it
    //   points to so the vector data will still be valid).
    must_drop: bool,
    phantom: PhantomData<&'a T>,
}

impl<T> Drop for View<'_, T> {
    fn drop(&mut self) {
        if self.must_drop {
            unsafe { ManuallyDrop::drop(&mut self.inner) };
        }
    }
}

impl<T> Deref for View<'_, T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T: Debug> Debug for View<'_, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        <T as Debug>::fmt(self, f)
    }
}

impl<'a, T> View<'a, T> {
    pub(crate) fn new(t: T, must_drop: bool) -> Self {
        Self {
            inner: ManuallyDrop::new(t),
            must_drop,
            phantom: PhantomData,
        }
    }
}

impl<'a, T: FFI + 'a> View<'a, T> {
    pub(crate) fn from_ptr(ptr: *const T::Sys, must_drop: bool) -> Self {
        // SAFETY: View only allows read-only access so the cast is OK.
        Self::new(FFI::wrap(ptr as *mut _), must_drop)
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
        Self::from_ptr(ptr, true)
    }
}

/// A mutable view to `T`.
pub struct ViewMut<'a, T> {
    inner: ManuallyDrop<T>,
    must_drop: bool,
    phantom: PhantomData<&'a ()>,
}

impl<T> Drop for ViewMut<'_, T> {
    fn drop(&mut self) {
        if self.must_drop {
            unsafe { ManuallyDrop::drop(&mut self.inner) };
        }
    }
}

impl<T> Deref for ViewMut<'_, T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T> DerefMut for ViewMut<'_, T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}

impl<T: FFI + Debug> Debug for ViewMut<'_, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        <T as Debug>::fmt(self, f)
    }
}

impl<'a, T> ViewMut<'a, T> {
    pub(crate) fn new(t: T, must_drop: bool) -> Self {
        Self {
            inner: ManuallyDrop::new(t),
            must_drop,
            phantom: PhantomData,
        }
    }
}

impl<'a, T: FFI + 'a> ViewMut<'a, T> {
    pub(crate) fn from_ptr(ptr: *mut T::Sys, must_drop: bool) -> Self {
        Self::new(FFI::wrap(ptr), must_drop)
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
        Self::from_ptr(ptr, true)
    }
}
