//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#![macro_use]

#[doc(hidden)]
macro_rules! ffi_wrap {
    ($name:tt) => {
        unsafe { $crate::ffi::FFI::wrap(sys::$name as *mut _) }
    };
}

#[doc(hidden)]
macro_rules! wrap_callback {
    ($f:expr, $F:ident $(+ $lt:lifetime)?) => {{
        unsafe extern "C" fn trampoline<$($lt,)? F: Fn(f64) -> f64 $( + $lt)?>(
            x: f64,
            params: *mut std::os::raw::c_void,
        ) -> f64 {
            let f: &F = unsafe { &*(params as *const F) };
            let x = f(x);
            x
        }

        sys::gsl_function_struct {
            function: Some(trampoline::<$F>),
            params: &$f as *const _ as *mut _,
        }
    }};
}

#[doc(hidden)]
macro_rules! ffi_wrapper {
    ($(#[$doc: meta])* $name:ident
        $(<$($p:tt $(: $tr0:tt $(+ ?$tr1:tt)?)?),*>)?,
        *mut $ty:ty, $drop:ident
        $(;$extra_id:ident: $extra_ty:ty => $extra_expr:expr;)*
    ) => {
        ffi_wrapper!($(#[$doc])* $name
            $(<$($p $(: $tr0 $(+ ?$tr1)?)?),*>)?,
            *mut $ty
            $(;$extra_id: $extra_ty => $extra_expr;)*);

        impl$(<$($p $(: $tr0 $(+ ?$tr1)?)?),*>)? Drop for $name$(<$($p),*>)? {
            fn drop(&mut self) {
                unsafe { sys::$drop(self.inner) };
                self.inner = std::ptr::null_mut();
            }
        }
    };
    ($(#[$doc:meta])* $name:ident
        $(<$($p:tt $(: $tr0:tt $(+ ?$tr1:tt)?)?),*>)?, *mut $ty:ty
        $(;$extra_id:ident: $extra_ty:ty => $extra_expr:expr;)*
    ) => {
        $(#[$doc])*
        pub struct $name$(<$($p $(: $tr0 $(+ ?$tr1)?)?),*>)? {
            inner: *mut $ty,
            $($extra_id: $extra_ty,)*
        }

        unsafe impl$(<$($p $(: $tr0 $(+ ?$tr1)?)?),*>)? FFI
        for $name$(<$($p),*>)? {
            type Sys = $ty;

            fn wrap(inner: *mut $ty) -> Self {
                Self { inner $(, $extra_id: $extra_expr)* }
            }

            #[inline]
            fn unwrap_shared(&self) -> *const $ty {
                self.inner as *const _
            }

            #[inline]
            fn unwrap_unique(&mut self) -> *mut $ty {
                self.inner
            }
        }
    };
    ($(#[$doc:meta])* $name:ident
        $(<$($p:tt $(: $tr0:tt $(+ ?$tr1:tt)?)?),*>)?, *const $ty:ty
        $(;$extra_id:ident: $extra_ty:ty => $extra_expr:expr;)*
    ) => {
        $(#[$doc])*
        #[derive(Clone, Copy)]
        pub struct $name$(<$($p $(: $tr0 $(+ ?$tr1)?)?),*>)? {
            inner: *const $ty,
            $($extra_id: $extra_ty,)*
        }

        unsafe impl$(<$($p$(: $tr0 $(+ ?$tr1)?)?),*>)? FFI
        for $name$(<$($p),*>)? {
            type Sys = $ty;

            fn wrap(inner: *mut $ty) -> Self {
                Self { inner $(, $extra_id: $extra_expr)* }
            }

            #[inline]
            fn unwrap_shared(&self) -> *const $ty {
                self.inner
            }

            fn unwrap_unique(&mut self) -> *mut $ty {
                unimplemented!()
            }
        }
    };
}

#[doc(hidden)]
macro_rules! map_name {
    ($fn_name:path, $map:expr, $name: ident, $T:ty) => {{
        use std::{collections::HashMap, ffi::CStr, sync::LazyLock};

        static MAP: LazyLock<HashMap<&CStr, $T>> = LazyLock::new(|| {
            let mut m = HashMap::new();
            for (n, t) in $map {
                m.insert(n, t);
            }
            m
        });
        if $name.is_null() {
            panic!("{}::name: null pointer", stringify!($fn_name));
        }
        let name: &CStr = unsafe { CStr::from_ptr($name) };
        *MAP.get(name).expect(&format!(
            "{}::name: {:?} unknown",
            stringify!($fn_name),
            name
        ))
    }};
}
