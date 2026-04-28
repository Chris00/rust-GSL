//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// # Safety
///
/// `Drop` must be implemented for `Self` in order to clear the
/// resources held by `*mut Self::Sys`.
#[allow(clippy::upper_case_acronyms)]
pub unsafe trait FFI {
    type Sys;

    /// Wrap the GSL pointer.
    fn wrap(r: *mut Self::Sys) -> Self;

    fn unwrap_shared(&self) -> *const Self::Sys;
    fn unwrap_unique(&mut self) -> *mut Self::Sys;
}
