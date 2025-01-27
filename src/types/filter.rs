//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::ffi::FFI;
use crate::{Error, FilterEnd, FilterScale, VectorF64, VectorI32};

ffi_wrapper!(
    FilterGaussianWorkspace,
    *mut sys::gsl_filter_gaussian_workspace,
    gsl_filter_gaussian_free
);

impl FilterGaussianWorkspace {
    #[doc(alias = "gsl_filter_gaussian_alloc")]
    pub fn new(K: usize) -> Option<Self> {
        let s = unsafe { sys::gsl_filter_gaussian_alloc(K) };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    /// This function applies a Gaussian filter parameterized by `alpha` to the input vector `x`,
    /// storing the output in `y`. The derivative order is specified by `order`, with `0`
    /// corresponding to a Gaussian, `1` corresponding to a first derivative Gaussian, and so on.
    /// The parameter `endtype` specifies how the signal end points are handled. It is allowed for
    /// `x` = `y` for an in-place filter.
    #[doc(alias = "gsl_filter_gaussian")]
    pub fn gaussian(
        &mut self,
        endtype: FilterEnd,
        alpha: f64,
        order: usize,
        x: &VectorF64,
        y: &mut VectorF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_filter_gaussian(
                endtype.into(),
                alpha,
                order,
                x.unwrap_shared(),
                y.unwrap_unique(),
                self.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
}

ffi_wrapper!(
    FilterMedianWorkspace,
    *mut sys::gsl_filter_median_workspace,
    gsl_filter_median_free
);

impl FilterMedianWorkspace {
    #[doc(alias = "gsl_filter_median_alloc")]
    pub fn new(K: usize) -> Option<Self> {
        let s = unsafe { sys::gsl_filter_median_alloc(K) };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    #[doc(alias = "gsl_filter_median")]
    pub fn median(
        &mut self,
        endtype: FilterEnd,
        x: &VectorF64,
        y: &mut VectorF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_filter_median(
                endtype.into(),
                x.unwrap_shared(),
                y.unwrap_unique(),
                self.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
}

ffi_wrapper!(
    FilterRMedianWorkspace,
    *mut sys::gsl_filter_rmedian_workspace,
    gsl_filter_rmedian_free
);

impl FilterRMedianWorkspace {
    #[doc(alias = "gsl_filter_rmedian_alloc")]
    pub fn new(K: usize) -> Option<Self> {
        let s = unsafe { sys::gsl_filter_rmedian_alloc(K) };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    #[doc(alias = "gsl_filter_rmedian")]
    pub fn rmedian(
        &mut self,
        endtype: FilterEnd,
        x: &VectorF64,
        y: &mut VectorF64,
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_filter_rmedian(
                endtype.into(),
                x.unwrap_shared(),
                y.unwrap_unique(),
                self.unwrap_unique(),
            )
        };
        Error::handle(ret, ())
    }
}

ffi_wrapper!(
    FilterImpulseWorkspace,
    *mut sys::gsl_filter_impulse_workspace,
    gsl_filter_impulse_free
);

impl FilterImpulseWorkspace {
    #[doc(alias = "gsl_filter_impulse_alloc")]
    pub fn new(K: usize) -> Option<Self> {
        let s = unsafe { sys::gsl_filter_impulse_alloc(K) };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    /// Returns `noutlier`.
    #[doc(alias = "gsl_filter_impulse")]
    pub fn impulse(
        &mut self,
        endtype: FilterEnd,
        scale_type: FilterScale,
        t: f64,
        x: &VectorF64,
        y: &mut VectorF64,
        xmedian: &mut VectorF64,
        xsigma: &mut VectorF64,
        ioutlier: &mut VectorI32,
    ) -> Result<usize, Error> {
        let mut noutlier = 0;
        let ret = unsafe {
            sys::gsl_filter_impulse(
                endtype.into(),
                scale_type.into(),
                t,
                x.unwrap_shared(),
                y.unwrap_unique(),
                xmedian.unwrap_unique(),
                xsigma.unwrap_unique(),
                &mut noutlier,
                ioutlier.unwrap_unique(),
                self.unwrap_unique(),
            )
        };
        Error::handle(ret, noutlier)
    }
}
