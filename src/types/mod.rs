//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! GLS types (reexported into the root `rgsl`).

pub use self::basis_spline::BSpLineWorkspace;

pub use self::chebyshev::ChebSeries;
pub use self::combination::Combination;
#[allow(deprecated)]
#[cfg(feature = "complex")]
pub use self::complex::{ComplexF32, ComplexF64, ComplexOps};
pub use self::discrete_hankel::DiscreteHankel;
pub use self::eigen_symmetric_workspace::{
    EigenGenHermVWorkspace, EigenGenHermWorkspace, EigenGenSymmVWorkspace, EigenGenSymmWorkspace,
    EigenGenVWorkspace, EigenGenWorkspace, EigenHermitianVWorkspace, EigenHermitianWorkspace,
    EigenNonSymmetricVWorkspace, EigenNonSymmetricWorkspace, EigenSymmetricVWorkspace,
    EigenSymmetricWorkspace,
};
#[cfg(feature = "complex")]
pub use self::fast_fourier_transforms::{
    FftComplexF32WaveTable, FftComplexF32Workspace, FftComplexF64WaveTable, FftComplexF64Workspace,
};
#[cfg(feature = "v2_5")]
#[cfg_attr(docsrs, doc(cfg(feature = "v2_5")))]
pub use self::filter::{
    FilterGaussianWorkspace, FilterImpulseWorkspace, FilterMedianWorkspace, FilterRMedianWorkspace,
};
pub use self::histograms::{Histogram, Histogram2D, Histogram2DPdf, HistogramPdf};
pub use self::interpolation::{Interp, InterpAccel, InterpType, Spline};
pub use self::matrix::{
    MatrixF32, MatrixF32View, MatrixF64, MatrixF64View, MatrixI32, MatrixI32View, MatrixU32,
    MatrixU32View,
};
#[cfg(feature = "complex")]
pub use self::matrix_complex::{MatrixComplexF32, MatrixComplexF64};
pub use self::minimizer::{Minimizer, MinimizerType};
pub use self::monte_carlo::{
    MiserMonteCarlo, MiserParams, PlainMonteCarlo, VegasMonteCarlo, VegasParams,
};
pub use self::multifit_linear::MultifitLinearWorkspace;
pub use self::multifit_solver::{
    MultiFitFdfSolver, MultiFitFdfSolverType, MultiFitFunction, MultiFitFunctionFdf,
};
#[cfg(feature = "v2_1")]
#[cfg_attr(docsrs, doc(cfg(feature = "v2_1")))]
pub use self::multilarge_linear::{MultilargeLinearType, MultilargeLinearWorkspace};
pub use self::multiset::MultiSet;
pub use self::n_tuples::{ReadNTuples, WriteNTuples};
pub use self::permutation::Permutation;
pub use self::polynomial::PolyComplexWorkspace;
pub use self::qrng::{QRng, QRngType};
pub use self::ran_discrete::RanDiscrete;
pub use self::result::{Result, ResultE10};
pub use self::rng::{Rng, RngType};
pub use self::roots::{RootFSolver, RootFSolverType, RootFdfSolver, RootFdfSolverType};
pub use self::rstat::{RStatQuantileWorkspace, RStatWorkspace};
pub use self::series_acceleration::{LevinUTruncWorkspace, LevinUWorkspace};
pub use self::siman::{SimAnnealing, SimAnnealingParams};
pub use self::vector::{
    VectorF32, VectorF32View, VectorF64, VectorF64View, VectorI32, VectorI32View, VectorU32,
    VectorU32View,
};
#[cfg(feature = "complex")]
pub use self::vector_complex::{
    VectorComplexF32, VectorComplexF32View, VectorComplexF64, VectorComplexF64View,
};
pub use self::wavelet_transforms::{Wavelet, WaveletType, WaveletWorkspace};

pub mod basis_spline;
pub mod chebyshev;
pub mod combination;
pub mod complex;
pub mod discrete_hankel;
pub mod eigen_symmetric_workspace;
pub mod fast_fourier_transforms;
#[cfg(feature = "v2_5")]
#[cfg_attr(docsrs, doc(cfg(feature = "v2_5")))]
pub mod filter;
pub mod histograms;
pub mod interpolation;
pub mod matrix;
pub mod matrix_complex;
pub mod minimizer;
pub mod monte_carlo;
pub mod multifit_linear;
pub mod multifit_solver;
#[cfg(feature = "v2_1")]
#[cfg_attr(docsrs, doc(cfg(feature = "v2_1")))]
pub mod multilarge_linear;
pub mod multimin;
pub mod multiroot;
pub mod multiset;
pub mod n_tuples;
pub mod permutation;
pub mod polynomial;
pub mod qrng;
pub mod ran_discrete;
pub mod result;
pub mod rng;
pub mod roots;
pub mod rstat;
pub mod series_acceleration;
pub mod siman;
pub mod vector;
pub mod vector_complex;
pub mod wavelet_transforms;
