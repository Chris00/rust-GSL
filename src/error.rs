//! Error handling.

use ctor::ctor;
use std::os::raw::c_int;

/// GSL errors.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Copy)]
pub enum Error {
    Failure,
    /// iteration has not converged
    Continue,
    /// input domain error, e.g sqrt(-1)
    Domain,
    /// output range error, e.g. exp(1e100)
    Range,
    /// invalid pointer
    Fault,
    /// invalid argument supplied by user
    Invalid,
    /// generic failure
    Failed,
    /// factorization failed
    Factorization,
    /// sanity check failed - shouldn't happen
    Sanity,
    /// malloc failed
    NoMemory,
    /// Problem with user-supplied function
    BadFunction,
    /// iterative process is out of control
    RunAway,
    /// exceeded max number of iterations
    MaxIteration,
    /// tried to divide by zero
    ZeroDiv,
    /// user specified an invalid tolerance
    BadTolerance,
    /// failed to reach the specified tolerance
    Tolerance,
    /// underflow
    UnderFlow,
    /// overflow
    OverFlow,
    /// loss of accuracy
    Loss,
    /// failed because of roundoff error
    Round,
    /// matrix, vector lengths are not conformant
    BadLength,
    /// matrix not square
    NotSquare,
    /// apparent singularity detected
    Singularity,
    /// integral or series is divergent
    Diverge,
    /// requested feature is not supported by the hardware
    Unsupported,
    /// requested feature not (yet) implemented
    Unimplemented,
    /// cache limit exceeded
    Cache,
    /// table limit exceeded
    Table,
    /// iteration is not making progress towards solution
    NoProgress,
    /// jacobian evaluations are not improving the solution
    NoProgressJacobian,
    /// cannot reach the specified tolerance in F
    ToleranceF,
    /// cannot reach the specified tolerance in X
    ToleranceX,
    /// cannot reach the specified tolerance in gradient
    ToleranceG,
    /// cannot reach the specified tolerance in gradient
    #[allow(clippy::upper_case_acronyms)]
    EOF,
    /// Unknown value.
    Unknown(i32),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use Error::*;
        let err = match self {
            Failure => "Failure",
            Continue => "The iteration has not converged yet",
            Domain => "Input domain error",
            Range => "Output range error",
            Fault => "Invalid pointer",
            Invalid => "Invalid argument supplied by user",
            Failed => "generic failure",
            Factorization => "Factorization failed",
            Sanity => "Sanity check failed - shouldn't happen",
            NoMemory => "Malloc failed",
            BadFunction => "Problem with user-supplied function",
            RunAway => "Iterative process is out of control",
            MaxIteration => "Exceeded max number of iterations",
            ZeroDiv => "Tried to divide by zero",
            BadTolerance => {
                "Specified tolerance is invalid or \
                            theoretically unattainable"
            }
            Tolerance => "Failed to reach the specified tolerance",
            UnderFlow => "Underflow",
            OverFlow => "Overflow",
            Loss => "Loss of accuracy",
            Round => "Roundoff error",
            BadLength => "Matrix/vector sizes are not conformant",
            NotSquare => "Matrix not square",
            Singularity => {
                "Singularity or extremely bad function \
                           behavior detected"
            }
            Diverge => "Integral or series is divergent",
            Unsupported => {
                "The required feature is not supported by \
                           this hardware platform"
            }
            Unimplemented => "The requested feature is not (yet) implemented",
            Cache => "Cache limit exceeded",
            Table => "Table limit exceeded",
            NoProgress => "Iteration is not making progress towards solution",
            NoProgressJacobian => {
                "Jacobian evaluations are not improving \
                                   the solution"
            }
            ToleranceF => "Cannot reach the specified tolerance in F",
            ToleranceX => "Cannot reach the specified tolerance in X",
            ToleranceG => "Cannot reach the specified tolerance in gradient",
            EOF => "End of file",
            Unknown(_) => "Unknown error",
        };
        write!(f, "{}", err)
    }
}

impl std::error::Error for Error {}

impl Error {
    pub(crate) fn handle<T>(v: c_int, x: T) -> Result<T, Error> {
        match v {
            sys::GSL_SUCCESS => Ok(x),
            sys::GSL_FAILURE => Err(Self::Failure),
            sys::GSL_CONTINUE => Err(Self::Continue),
            sys::GSL_EDOM => Err(Self::Domain),
            sys::GSL_ERANGE => Err(Self::Range),
            sys::GSL_EFAULT => Err(Self::Fault),
            sys::GSL_EINVAL => Err(Self::Invalid),
            sys::GSL_EFAILED => Err(Self::Failed),
            sys::GSL_EFACTOR => Err(Self::Factorization),
            sys::GSL_ESANITY => Err(Self::Sanity),
            sys::GSL_ENOMEM => Err(Self::NoMemory),
            sys::GSL_EBADFUNC => Err(Self::BadFunction),
            sys::GSL_ERUNAWAY => Err(Self::RunAway),
            sys::GSL_EMAXITER => Err(Self::MaxIteration),
            sys::GSL_EZERODIV => Err(Self::ZeroDiv),
            sys::GSL_EBADTOL => Err(Self::BadTolerance),
            sys::GSL_ETOL => Err(Self::Tolerance),
            sys::GSL_EUNDRFLW => Err(Self::UnderFlow),
            sys::GSL_EOVRFLW => Err(Self::OverFlow),
            sys::GSL_ELOSS => Err(Self::Loss),
            sys::GSL_EROUND => Err(Self::Round),
            sys::GSL_EBADLEN => Err(Self::BadLength),
            sys::GSL_ENOTSQR => Err(Self::NotSquare),
            sys::GSL_ESING => Err(Self::Singularity),
            sys::GSL_EDIVERGE => Err(Self::Diverge),
            sys::GSL_EUNSUP => Err(Self::Unsupported),
            sys::GSL_EUNIMPL => Err(Self::Unimplemented),
            sys::GSL_ECACHE => Err(Self::Cache),
            sys::GSL_ETABLE => Err(Self::Table),
            sys::GSL_ENOPROG => Err(Self::NoProgress),
            sys::GSL_ENOPROGJ => Err(Self::NoProgressJacobian),
            sys::GSL_ETOLF => Err(Self::ToleranceF),
            sys::GSL_ETOLX => Err(Self::ToleranceX),
            sys::GSL_ETOLG => Err(Self::ToleranceG),
            sys::GSL_EOF => Err(Self::EOF),
            x => Err(Self::Unknown(x)),
        }
    }

    pub(crate) fn to_c(x: std::result::Result<(), Error>) -> c_int {
        match x {
            Ok(()) => sys::GSL_SUCCESS,
            Err(Error::Failure) => sys::GSL_FAILURE,
            Err(Error::Continue) => sys::GSL_CONTINUE,
            Err(Error::Domain) => sys::GSL_EDOM,
            Err(Error::Range) => sys::GSL_ERANGE,
            Err(Error::Fault) => sys::GSL_EFAULT,
            Err(Error::Invalid) => sys::GSL_EINVAL,
            Err(Error::Failed) => sys::GSL_EFAILED,
            Err(Error::Factorization) => sys::GSL_EFACTOR,
            Err(Error::Sanity) => sys::GSL_ESANITY,
            Err(Error::NoMemory) => sys::GSL_ENOMEM,
            Err(Error::BadFunction) => sys::GSL_EBADFUNC,
            Err(Error::RunAway) => sys::GSL_ERUNAWAY,
            Err(Error::MaxIteration) => sys::GSL_EMAXITER,
            Err(Error::ZeroDiv) => sys::GSL_EZERODIV,
            Err(Error::BadTolerance) => sys::GSL_EBADTOL,
            Err(Error::Tolerance) => sys::GSL_ETOL,
            Err(Error::UnderFlow) => sys::GSL_EUNDRFLW,
            Err(Error::OverFlow) => sys::GSL_EOVRFLW,
            Err(Error::Loss) => sys::GSL_ELOSS,
            Err(Error::Round) => sys::GSL_EROUND,
            Err(Error::BadLength) => sys::GSL_EBADLEN,
            Err(Error::NotSquare) => sys::GSL_ENOTSQR,
            Err(Error::Singularity) => sys::GSL_ESING,
            Err(Error::Diverge) => sys::GSL_EDIVERGE,
            Err(Error::Unsupported) => sys::GSL_EUNSUP,
            Err(Error::Unimplemented) => sys::GSL_EUNIMPL,
            Err(Error::Cache) => sys::GSL_ECACHE,
            Err(Error::Table) => sys::GSL_ETABLE,
            Err(Error::NoProgress) => sys::GSL_ENOPROG,
            Err(Error::NoProgressJacobian) => sys::GSL_ENOPROGJ,
            Err(Error::ToleranceF) => sys::GSL_ETOLF,
            Err(Error::ToleranceX) => sys::GSL_ETOLX,
            Err(Error::ToleranceG) => sys::GSL_ETOLG,
            Err(Error::EOF) => sys::GSL_EOF,
            Err(Error::Unknown(x)) => x,
        }
    }
}

// The GSL manual clearly says that the purpose of setting an error
// handler "is to provide a function where a breakpoint can be set
// that will catch library errors when running under the debugger.  It
// is not intended for use in production programs, which should handle
// any errors using the error return codes".  In Rust, errors cannot
// be ignored, so we deactivate the error handler.
#[ctor]
fn disable_error_handler() {
    unsafe {
        sys::gsl_set_error_handler_off();
    }
}
