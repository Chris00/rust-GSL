//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#[deprecated(since = "8.0.0", note = "Use rgsl::sf::Prec instead")]
pub type Mode = crate::sf::Prec;

#[deprecated(since = "8.0.0", note = "Use rgsl::Error instead")]
pub type Value = crate::Error;

#[deprecated(since = "8.0.0", note = "Use rgsl::eigen::Sort instead")]
pub type EigenSort = crate::eigen::Sort;

#[deprecated(since = "8.0.0", note = "Use rgsl::fft::Dir instead")]
pub type FftDirection = crate::fft::Dir;

/// Used by [`VegasParams`][crate::VegasParams].
///
/// This determines whether vegas will use importance sampling or
/// stratified sampling, or whether it can pick on its own.  In low
/// dimensions vegas uses strict stratified sampling (more precisely,
/// stratified sampling is chosen if there are fewer than 2 bins per
/// box).
#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum VegasMode {
    /// Importance sampling: allocate more sample points where the
    /// integrand is larger.
    Importance,
    /// Exclusively use importance sampling without any stratification.
    ImportanceOnly,
    /// Stratified sampling: divides the integration region into
    /// sub-regions and sample each sub-region separately.
    Stratified,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<std::os::raw::c_int> for VegasMode {
    fn into(self) -> std::os::raw::c_int {
        match self {
            Self::Importance => sys::GSL_VEGAS_MODE_IMPORTANCE,
            Self::ImportanceOnly => sys::GSL_VEGAS_MODE_IMPORTANCE_ONLY,
            Self::Stratified => sys::GSL_VEGAS_MODE_STRATIFIED,
        }
    }
}

#[doc(hidden)]
impl From<std::os::raw::c_int> for VegasMode {
    fn from(v: std::os::raw::c_int) -> VegasMode {
        match v {
            sys::GSL_VEGAS_MODE_IMPORTANCE => Self::Importance,
            sys::GSL_VEGAS_MODE_IMPORTANCE_ONLY => Self::ImportanceOnly,
            sys::GSL_VEGAS_MODE_STRATIFIED => Self::Stratified,
            _ => panic!("Unknown VegasMode value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum WaveletDirection {
    Forward,
    Backward,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::gsl_wavelet_direction> for WaveletDirection {
    fn into(self) -> sys::gsl_wavelet_direction {
        match self {
            Self::Forward => sys::gsl_wavelet_direction_gsl_wavelet_forward,
            Self::Backward => sys::gsl_wavelet_direction_gsl_wavelet_backward,
        }
    }
}

#[doc(hidden)]
impl From<sys::gsl_wavelet_direction> for WaveletDirection {
    fn from(v: sys::gsl_wavelet_direction) -> WaveletDirection {
        match v {
            sys::gsl_wavelet_direction_gsl_wavelet_forward => Self::Forward,
            sys::gsl_wavelet_direction_gsl_wavelet_backward => Self::Backward,
            _ => panic!("Unknown WaveletDirection value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum SfLegendreNorm {
    Schmidt,
    SphericalHarmonic,
    Full,
    None,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::gsl_sf_legendre_t> for SfLegendreNorm {
    fn into(self) -> sys::gsl_sf_legendre_t {
        match self {
            SfLegendreNorm::Schmidt => sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_SCHMIDT,
            SfLegendreNorm::SphericalHarmonic => sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_SPHARM,
            SfLegendreNorm::Full => sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_FULL,
            SfLegendreNorm::None => sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_NONE,
        }
    }
}

#[doc(hidden)]
impl From<sys::gsl_sf_legendre_t> for SfLegendreNorm {
    fn from(v: sys::gsl_sf_legendre_t) -> SfLegendreNorm {
        match v {
            sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_SCHMIDT => SfLegendreNorm::Schmidt,
            sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_SPHARM => SfLegendreNorm::SphericalHarmonic,
            sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_FULL => SfLegendreNorm::Full,
            sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_NONE => SfLegendreNorm::None,
            _ => panic!("Unknown SfLegendreNorm value"),
        }
    }
}

#[cfg(feature = "v2_5")]
#[cfg_attr(docsrs, doc(cfg(feature = "v2_5")))]
#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum FilterEnd {
    PadZero,
    PadValue,
    Truncate,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
#[cfg(feature = "v2_5")]
impl Into<sys::gsl_filter_end_t> for FilterEnd {
    fn into(self) -> sys::gsl_filter_end_t {
        match self {
            Self::PadZero => sys::gsl_filter_end_t_GSL_FILTER_END_PADZERO,
            Self::PadValue => sys::gsl_filter_end_t_GSL_FILTER_END_PADVALUE,
            Self::Truncate => sys::gsl_filter_end_t_GSL_FILTER_END_TRUNCATE,
        }
    }
}

#[doc(hidden)]
#[cfg(feature = "v2_5")]
impl From<sys::gsl_filter_end_t> for FilterEnd {
    fn from(v: sys::gsl_filter_end_t) -> FilterEnd {
        match v {
            sys::gsl_filter_end_t_GSL_FILTER_END_PADZERO => Self::PadZero,
            sys::gsl_filter_end_t_GSL_FILTER_END_PADVALUE => Self::PadValue,
            sys::gsl_filter_end_t_GSL_FILTER_END_TRUNCATE => Self::Truncate,
            _ => panic!("Unknown FilterEnd value"),
        }
    }
}

#[cfg(feature = "v2_5")]
#[cfg_attr(docsrs, doc(cfg(feature = "v2_5")))]
#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum FilterScale {
    MedianAbsoluteDeviation,
    InterQuartileRange,
    SN,
    QN,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
#[cfg(feature = "v2_5")]
impl Into<sys::gsl_filter_scale_t> for FilterScale {
    fn into(self) -> sys::gsl_filter_scale_t {
        match self {
            Self::MedianAbsoluteDeviation => sys::gsl_filter_scale_t_GSL_FILTER_SCALE_MAD,
            Self::InterQuartileRange => sys::gsl_filter_scale_t_GSL_FILTER_SCALE_IQR,
            Self::SN => sys::gsl_filter_scale_t_GSL_FILTER_SCALE_SN,
            Self::QN => sys::gsl_filter_scale_t_GSL_FILTER_SCALE_QN,
        }
    }
}

#[doc(hidden)]
#[cfg(feature = "v2_5")]
impl From<sys::gsl_filter_scale_t> for FilterScale {
    fn from(v: sys::gsl_filter_scale_t) -> FilterScale {
        match v {
            sys::gsl_filter_scale_t_GSL_FILTER_SCALE_MAD => Self::MedianAbsoluteDeviation,
            sys::gsl_filter_scale_t_GSL_FILTER_SCALE_IQR => Self::InterQuartileRange,
            sys::gsl_filter_scale_t_GSL_FILTER_SCALE_SN => Self::SN,
            sys::gsl_filter_scale_t_GSL_FILTER_SCALE_QN => Self::QN,
            _ => panic!("Unknown FilterScale value"),
        }
    }
}
