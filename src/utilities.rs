//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use std::ffi::CString;
use std::io;
use std::os::raw::c_char;
use std::path::Path;

use libc::{FILE, fclose, fopen};

/// # Safety
/// The caller must ensure that `f` lives for as long at the return
/// value is used.
pub(crate) unsafe fn wrap_callback<F: FnMut(f64) -> f64>(f: &mut F) -> sys::gsl_function_struct {
    unsafe extern "C" fn trampoline<F: FnMut(f64) -> f64>(
        x: f64,
        params: *mut std::os::raw::c_void,
    ) -> f64 {
        let f = unsafe { &mut *params.cast::<F>() };
        f(x)
    }

    sys::gsl_function_struct {
        function: Some(trampoline::<F>),
        params: &mut *f as *const _ as *mut _,
    }
}

/// Return a GSL callback struct that is allocated on the heap.  This
/// is necessary in order to pass a pointer to
/// `sys::gsl_function_struct` as an argument, as the address must
/// remain fixed across several calls during which the object holding
/// the returned value of this function (so it is not dropped) may be
/// moved.
///
/// # Safety
/// The caller must ensure that `f` lives for as long at the return
/// value is used.
pub(crate) unsafe fn box_callback<F: FnMut(f64) -> f64>(
    f: &mut Box<F>,
) -> Box<sys::gsl_function_struct> {
    Box::new(unsafe { wrap_callback(&mut *f) })
}

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mode {
    Write,
    Read,
}

/// A wrapper to handle I/O operations between GSL and rust
#[allow(clippy::upper_case_acronyms)]
pub struct IOStream {
    inner: *mut FILE,
    mode: Mode,
}

impl IOStream {
    /// Open a file in write mode.
    pub fn fwrite_handle<P: AsRef<Path>>(file: &P) -> io::Result<IOStream> {
        let path = CString::new(file.as_ref().to_str().unwrap()).unwrap();
        let ptr = unsafe { fopen(path.as_ptr(), c"w".as_ptr() as *const c_char) };
        if ptr.is_null() {
            return Err(io::Error::other("Failed to open file..."));
        }
        Ok(IOStream {
            inner: ptr,
            mode: Mode::Write,
        })
    }

    pub fn write_mode(&self) -> bool {
        self.mode == Mode::Write
    }

    #[doc(hidden)]
    pub fn as_raw(&mut self) -> *mut FILE {
        self.inner
    }
}

impl Drop for IOStream {
    fn drop(&mut self) {
        unsafe {
            fclose(self.inner);
            self.inner = std::ptr::null_mut();
        }
    }
}
