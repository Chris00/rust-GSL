//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Nonlinear Least-Squares Fitting

This module contains functions for multidimensional nonlinear
least-squares fitting.  There are generally two classes of algorithms
for solving nonlinear least squares problems, which fall under line
search methods and trust region methods.  GSL currently implements
only trust region methods and provides the user with full access to
intermediate steps of the iteration.  The user also has the ability to
tune a number of parameters which affect low-level aspects of the
algorithm which can help to accelerate convergence for the specific
problem at hand.  GSL provides two separate interfaces for nonlinear
least squares fitting.  The first is designed for small to moderate
sized problems, and the second is designed for very large problems,
which may or may not have significant sparse structure.

This module [`rgsl::multifit`][crate::multifit] contains the functions
for the multidimensional nonlinear fitting functions and related
declarations relating to the small to moderate sized systems.

The submodule [`rgsl::multifit::large`][crate::multifit::large],
aliased to [`rgsl::multilarge`][crate::multilarge] contains the
functions for the multidimensional nonlinear fitting functions and
related declarations relating to large systems.

## Overview

The problem of multidimensional nonlinear least-squares fitting
requires the minimization of the squared residuals of $n$ functions,
$f_i$, in $p$ parameters, $x_i$,

$$Φ(x) = \frac{1}{2} ‖f(x)‖²
= \frac{1}{2} ∑_{i=1}^n f_i (x_1,…, x_p)²$$

In trust region methods, the objective (or cost) function $Φ(x)$ is
approximated by a model function $mₖ(δ)$ in the vicinity of some
point $xₖ$.  The model function is often simply a second order Taylor
series expansion around the point $xₖ$, i.e.:

$$Φ(xₖ + δ) ≈ mₖ(\delta) = Φ(xₖ) + gₖ^T δ + \frac{1}{2} δ^T Bₖ δ$$

where $gₖ = ∇Φ(xₖ) = J^T f$ is the gradient vector at the point $xₖ$,
$Bₖ = ∇²Φ(xₖ)$ is the Hessian matrix at $xₖ$, or some approximation to
it, and $J$ is the $n$-by-$p$ Jacobian matrix

$$J_{ij} = ∂fᵢ / ∂xⱼ$$

In order to find the next step $δ$, we minimize the model function
$mₖ(δ)$, but search for solutions only within a region where we trust
that $mₖ(δ)$ is a good approximation to the objective function
$Φ(xₖ+δ)$.  In other words, we seek a solution of the trust region
subproblem ([`TRS`])

$$\min_{δ ∈ Rᵖ} mₖ(δ) = Φ(xₖ) + gₖ^T δ + \frac{1}{2} δ^T Bₖ δ,
\qquad\text{s.t.}\quad ‖Dₖ δ‖ ≤ Δₖ$$

where $Δₖ > 0$ is the trust region radius and $Dₖ$ is a scaling
matrix.  If $Dₖ = I$, then the trust region is a ball of radius $Δₖ$
centered at $xₖ$.  In some applications, the parameter vector $x$ may
have widely different scales.  For example, one parameter might be a
temperature on the order of $10³$ K, while another might be a length
on the order of $10^{-6}$ m.  In such cases, a spherical trust region
may not be the best choice, since if $Φ$ changes rapidly along
directions with one scale, and more slowly along directions with a
different scale, the model function $mₖ$ may be a poor approximation
to $Φ$ along the rapidly changing directions.  In such problems, it
may be best to use an elliptical trust region, by setting $Dₖ$ to a
diagonal matrix whose entries are designed so that the scaled step
$Dₖ δ$ has entries of approximately the same order of magnitude.

The trust region subproblem above normally amounts to solving a linear
least squares system (or multiple systems) for the step $δ$.  Once $δ$
is computed, it is checked whether or not it reduces the objective
function $Φ(x)$.  A useful statistic for this is to look at the ratio

$$ρₖ = \frac{Φ(xₖ) - Φ(xₖ + δₖ)}{mₖ(0) - mₖ(δₖ)}$$

where the numerator is the actual reduction of the objective function
due to the step $δₖ$, and the denominator is the predicted
reduction due to the model $mₖ$. If $ρₖ$ is negative, it means that
the step $δₖ$ increased the objective function and so it is
rejected.  If $ρₖ$ is positive, then we have found a step which
reduced the objective function and it is accepted.  Furthermore, if
$ρₖ$ is close to 1, then this indicates that the model function is a
good approximation to the objective function in the trust region, and
so on the next iteration the trust region is enlarged in order to take
more ambitious steps.  When a step is rejected, the trust region is
made smaller and the TRS is solved again.  An outline for the general
trust region method used by GSL can now be given.

**Trust Region Algorithm**

1. Initialize: given $x₀$, construct $m₀(δ)$, $D₀$ and $Δ₀ > 0$.

2. For $k = 0, 1, 2, …$

   a. If converged, then stop

   b. Solve TRS for trial step $δₖ$

   c. Evaluate trial step by computing $ρₖ$
      - if step is accepted, set $x_{k+1} = xₖ + δₖ$ and increase radius,
        $Δ_{k+1} = α Δₖ$
      - if step is rejected, set $x_{k+1} = xₖ$ and decrease radius,
        $Δ_{k+1} = Δₖ/β$; goto 2(b)

   d. Construct $m_{k+1}(δ)$ and $D_{k+1}$

GSL offers the user a number of different algorithms for solving the
trust region subproblem in 2(b), as well as different choices of
scaling matrices $Dₖ$ and different methods of updating the trust
region radius $Δₖ$.  Therefore, while reasonable default methods are
provided, the user has a lot of control to fine-tune the various steps
of the algorithm for their specific problem.

## Solving the Trust Region Subproblem (TRS)

Below we describe the methods available for solving the trust region
subproblem.  The methods available provide either exact or approximate
solutions to the trust region subproblem.  In all algorithms below,
the Hessian matrix $Bₖ$ is approximated as $Bₖ ≈ Jₖ^T Jₖ$, where $Jₖ =
J(xₖ)$.  In all methods, the solution of the TRS involves solving a
linear least squares system involving the Jacobian matrix.  For small
to moderate sized problems ([`MultiFitFSolver`] interface), this is
accomplished by factoring the full Jacobian matrix, which is provided
by the user, with the Cholesky, QR, or SVD decompositions.  For large
systems (`large::MultiFitFSolver` interface), the user has two
choices.  One is to solve the system iteratively, without needing to
store the full Jacobian matrix in memory.  With this method, the user
must provide a routine to calculate the matrix-vector products $J u$
or $J^T u$ for a given vector $u$.  This iterative method is
particularly useful for systems where the Jacobian has sparse
structure, since forming matrix-vector products can be done cheaply.
The second option for large systems involves forming the normal
equations matrix $J^T J$ and then factoring it using a Cholesky
decomposition.  The normal equations matrix is $p$-by-$p$, typically
much smaller than the full $n$-by-$p$ Jacobian, and can usually be
stored in memory even if the full Jacobian matrix cannot.  This option
is useful for large, dense systems, or if the iterative method has
difficulty converging.

### Levenberg-Marquardt

There is a theorem which states that if $δₖ$ is a solution to the
trust region subproblem given above, then there exists $μₖ ≥ 0$ such
that

$$\left( Bₖ + μₖ Dₖ^T Dₖ \right) δₖ = -gₖ$$

with $μₖ (Δₖ - ‖Dₖ δₖ‖) = 0$.  This forms the basis of the
Levenberg-Marquardt algorithm, which controls the trust region size by
adjusting the parameter $μₖ$ rather than the radius $Δₖ$ directly.
For each radius $Δₖ$, there is a unique parameter $μₖ$ which solves
the TRS, and they have an inverse relationship, so that large values
of $μₖ$ correspond to smaller trust regions, while small values of
$μₖ$ correspond to larger trust regions.

With the approximation $Bₖ ≈ Jₖ^T Jₖ$, on each iteration, in order to
calculate the step $δₖ$, the following linear least squares problem is
solved:

$$\begin{pmatrix} J_k \cr \sqrt{μₖ} Dₖ \end{pmatrix} \delta_k
= - \begin{pmatrix} fₖ \cr 0 \end{pmatrix}$$

If the step $δₖ$ is accepted, then $μₖ$ is decreased on the next
iteration in order to take a larger step, otherwise it is increased to
take a smaller step.  The Levenberg-Marquardt algorithm provides an
exact solution of the trust region subproblem, but typically has a
higher computational cost per iteration than the approximate methods
discussed below, since it may need to solve the least squares system
above several times for different values of $μₖ$.

### Levenberg-Marquardt with Geodesic Acceleration

This method applies a so-called geodesic acceleration correction to
the standard Levenberg-Marquardt step $δₖ$ (Transtrum et al, 2011).
By interpreting $δₖ$ as a first order step along a geodesic in the
model parameter space (i.e. a velocity $δₖ = vₖ$), the geodesic
acceleration $aₖ$ is a second order correction along the geodesic which
is determined by solving the linear least squares system

$$\begin{pmatrix} J_k \cr √{μₖ} Dₖ \end{pmatrix} aₖ
= - \begin{pmatrix} f_{vv}(xₖ) \cr 0 \end{pmatrix}$$

where $f_{vv}$ is the second directional derivative of the residual
vector in the velocity direction $v$, $f_{vv}(x) = Dᵥ² f = ∑_{αβ} v_α
v_β ∂_α ∂_β f(x)$, where $α$ and $β$ are summed over the $p$
parameters.  The new total step is then $δₖ′ = vₖ + \frac{1}{2}aₖ$.
The second order correction $aₖ$ can be calculated with a modest
additional cost, and has been shown to dramatically reduce the number
of iterations (and expensive Jacobian evaluations) required to reach
convergence on a variety of different problems.  In order to utilize
the geodesic acceleration, the user must supply a function which
provides the second directional derivative vector $f\_{vv}(x)$, or
alternatively the library can use a finite difference method to
estimate this vector with one additional function evaluation of
$f(x+hv)$ where $h$ is a tunable step size (see the
[`h_fvv`](Parameters::h_fvv) parameter description).

### Dogleg

This is Powell’s dogleg method, which finds an approximate solution to
the trust region subproblem, by restricting its search to a piecewise
linear “dogleg” path, composed of the origin, the Cauchy point which
represents the model minimizer along the steepest descent direction,
and the Gauss-Newton point, which is the overall minimizer of the
unconstrained model.  The Gauss-Newton step is calculated by solving

$$J_k δ_{gn} = -f_k$$

which is the main computational task for each iteration, but only
needs to be performed once per iteration.  If the Gauss-Newton point
is inside the trust region, it is selected as the step.  If it is
outside, the method then calculates the Cauchy point, which is located
along the gradient direction.  If the Cauchy point is also outside the
trust region, the method assumes that it is still far from the minimum
and so proceeds along the gradient direction, truncating the step at
the trust region boundary.  If the Cauchy point is inside the trust
region, with the Gauss-Newton point outside, the method uses a dogleg
step, which is a linear combination of the gradient direction and the
Gauss-Newton direction, stopping at the trust region boundary.

### Double Dogleg

This method is an improvement over the classical dogleg algorithm,
which attempts to include information about the Gauss-Newton step
while the iteration is still far from the minimum. When the Cauchy
point is inside the trust region and the Gauss-Newton point is
outside, the method computes a scaled Gauss-Newton point and then
takes a dogleg step between the Cauchy point and the scaled
Gauss-Newton point.  The scaling is calculated to ensure that the
 reduction in the model $mₖ$ is about the same as the reduction
provided by the Cauchy point.

### Two Dimensional Subspace

The dogleg methods restrict the search for the TRS solution to a 1D
curve defined by the Cauchy and Gauss-Newton points.  An improvement
to this is to search for a solution using the full two dimensional
subspace spanned by the Cauchy and Gauss-Newton directions.  The
dogleg path is of course inside this subspace, and so this method
solves the TRS at least as accurately as the dogleg methods.  Since
this method searches a larger subspace for a solution, it can converge
more quickly than dogleg on some problems.  Because the subspace is
only two dimensional, this method is very efficient and the main
computation per iteration is to determine the Gauss-Newton point.

### Steihaug-Toint Conjugate Gradient

One difficulty of the dogleg methods is calculating the Gauss-Newton
step when the Jacobian matrix is singular.  The Steihaug-Toint method
also computes a generalized dogleg step, but avoids solving for the
Gauss-Newton step directly, instead using an iterative conjugate
gradient algorithm.  This method performs well at points where the
Jacobian is singular, and is also suitable for large-scale problems
where factoring the Jacobian matrix could be prohibitively expensive.


## Weighted Nonlinear Least-Squares

Weighted nonlinear least-squares fitting minimizes the function

$$Φ(x) = \frac{1}{2} ‖f‖²_W
= \frac{1}{2} ∑_{i=1}^n wᵢ fᵢ (x₁,…, xₚ)²$$

where $W = \diag(w_1,w_2,...,w_n)$ is the weighting matrix, and
$‖f‖²_W = f^T W f$.  The weights $wᵢ$ are commonly defined as $wᵢ =
1/σᵢ²$, where $σᵢ$ is the error in the $i$-th measurement.  A simple
change of variables $\tilde{f} = W^{1/2} f$ yields $Φ(x) = \frac{1}{2}
‖\tilde{f}‖²$, which is in the same form as the unweighted case.  The
user can either perform this transform directly on their function
residuals and Jacobian, or use the `gsl_multifit_nlinear_winit()`
interface which automatically performs the correct scaling.  To
manually perform this transformation, the residuals and Jacobian
should be modified according to

$$\tilde{f}ᵢ = √{wᵢ} fᵢ = \frac{fᵢ}{σᵢ}$$
$$\tilde{J}_{ij} = √{wᵢ} \frac{∂fᵢ}{∂xⱼ}
= \frac{1}{σᵢ} \frac{∂fᵢ}{∂xⱼ}$$

For large systems, the user must perform their own weighting.

## Tunable Parameters

The user can tune nearly all aspects of the iteration at allocation
time.  See [`Parameters`] and [`large::Parameters`].
 */

use crate::{
    Error,
    ffi::FFI,
    matrix::{AsMatrix, MatF64},
    vector::VecF64,
};
use pastey::paste;
use std::ffi::c_void;

// `sqrt` is not const.
const SQRT_EPSILON: f64 = 1.4901161193847656e-08;

/// Method to solve the Trust Region Subproblem.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TRS {
    /// [Levenberg-Marquardt](crate::multifit#levenberg-marquardt)
    /// algorithm.
    LM,
    /// [Levenberg-Marquardt with Geodesic
    /// Acceleration](crate::multifit#levenberg-marquardt-with-geodesic-acceleration).
    LMaccel,
    /// [Dogleg](crate::multifit#dogleg) algorithm.
    Dogleg,
    /// [Double Dogleg](crate::multifit#double-dogleg) algorithm.
    DDogleg,
    /// [Two Dimensional
    /// Subspace](crate::multifit#two-dimensional-subspace) algorithm.
    Subspace2D,
}

impl TRS {
    fn to_c(self) -> *const sys::gsl_multifit_nlinear_trs {
        unsafe {
            match self {
                Self::LM => sys::gsl_multifit_nlinear_trs_lm,
                Self::LMaccel => sys::gsl_multifit_nlinear_trs_lmaccel,
                Self::Dogleg => sys::gsl_multifit_nlinear_trs_dogleg,
                Self::DDogleg => sys::gsl_multifit_nlinear_trs_ddogleg,
                Self::Subspace2D => sys::gsl_multifit_nlinear_trs_subspace2D,
            }
        }
    }
}

/// Determines the [diagonal scaling matrix
/// $D$](crate::multifit#overview).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scale {
    /// This damping strategy was suggested by Moré, and corresponds
    /// to $D^T D = \max(\diag(J^T J))$, in other words the maximum
    /// elements of $\diag(J^T J)$ encountered thus far in the
    /// iteration.  This choice of $D$ makes the problem
    /// scale-invariant, so that if the model parameters $xᵢ$ are each
    /// scaled by an arbitrary constant, $\tilde{x}ᵢ = aᵢ xᵢ$, then
    /// the sequence of iterates produced by the algorithm would be
    /// unchanged.  This method can work very well in cases where the
    /// model parameters have widely different scales (i.e. if some
    /// parameters are measured in nanometers, while others are
    /// measured in degrees Kelvin).  This strategy has been proven
    /// effective on a large class of problems and so it is the
    /// library default, but it may not be the best choice for all
    /// problems.
    More,
    /// This damping strategy was originally suggested by Levenberg,
    /// and corresponds to $D^T D = I$.  This method has also proven
    /// effective on a large class of problems, but is not
    /// scale-invariant.  However, some authors (e.g. Transtrum and
    /// Sethna 2012) argue that this choice is better for problems
    /// which are susceptible to parameter evaporation
    /// (i.e. parameters go to infinity)
    Levenberg,
    /// This damping strategy was suggested by Marquardt, and
    /// corresponds to $D^T D = \diag(J^T J)$.  This method is
    /// scale-invariant, but it is generally considered inferior to
    /// both the Levenberg and Moré strategies, though may work well
    /// on certain classes of problems.
    Marquardt,
}

impl Scale {
    fn to_c(self) -> *const sys::gsl_multifit_nlinear_scale {
        unsafe {
            match self {
                Self::More => sys::gsl_multifit_nlinear_scale_more,
                Self::Levenberg => sys::gsl_multifit_nlinear_scale_levenberg,
                Self::Marquardt => sys::gsl_multifit_nlinear_scale_marquardt,
            }
        }
    }

    fn to_c_large(self) -> *const sys::gsl_multilarge_nlinear_scale {
        unsafe {
            match self {
                Self::More => sys::gsl_multilarge_nlinear_scale_more,
                Self::Levenberg => sys::gsl_multilarge_nlinear_scale_levenberg,
                Self::Marquardt => sys::gsl_multilarge_nlinear_scale_marquardt,
            }
        }
    }
}

/// Solver for the trust region subproblem.
///
/// Solving the trust region subproblem on each iteration almost
/// always requires the solution of the following linear least squares
/// system
///
/// $$\begin{pmatrix} J \cr √μ D \end{pmatrix} δ
/// = - \begin{pmatrix} f \cr 0 \end{pmatrix}$$
///
/// How the system is solved and can be selected the choices of this
/// enum.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Solver {
    /// This method solves the system using a rank revealing QR
    /// decomposition of the Jacobian $J$.  This method will produce
    /// reliable solutions in cases where the Jacobian is rank
    /// deficient or near-singular but does require about twice as
    /// many operations as the `Cholesky` method.
    QR,
    /// This method solves the alternate normal equations problem
    ///
    /// $$\left( J^T J + μ D^T D \right) δ = -J^T f$$
    ///
    /// by using a Cholesky decomposition of the matrix $J^T J + μ D^T
    /// D$.  This method is faster than the QR approach, however it is
    /// susceptible to numerical instabilities if the Jacobian matrix
    /// is rank deficient or near-singular.  In these cases, an
    /// attempt is made to reduce the condition number of the matrix
    /// using Jacobi preconditioning, but for highly ill-conditioned
    /// problems the QR approach is better.  If it is known that the
    /// Jacobian matrix is well conditioned, this method is accurate
    /// and will perform faster than the QR approach.
    Cholesky,
    /// This method solves the alternate normal equations problem
    ///
    /// $$\left( J^T J + μ D^T D \right) δ = -J^T f$$
    ///
    /// by using a modified Cholesky decomposition of the matrix $J^T J
    /// + μ D^T D$.  This is more suitable for the dogleg methods
    /// where the parameter $μ = 0$, and the matrix $J^T J$ may be
    /// ill-conditioned or indefinite causing the standard Cholesky
    /// decomposition to fail.  This method is based on Level 2 BLAS
    /// and is thus slower than the standard Cholesky decomposition,
    /// which is based on Level 3 BLAS.
    MCholesky,
    /// This method solves the system using a singular value
    /// decomposition of the Jacobian $J$.  This method will produce
    /// the most reliable solutions for ill-conditioned Jacobians but
    /// is also the slowest solver method.
    SVD,
}

impl Solver {
    fn to_c(self) -> *const sys::gsl_multifit_nlinear_solver {
        unsafe {
            match self {
                Self::QR => sys::gsl_multifit_nlinear_solver_qr,
                Self::Cholesky => sys::gsl_multifit_nlinear_solver_cholesky,
                Self::MCholesky => sys::gsl_multifit_nlinear_solver_mcholesky,
                Self::SVD => sys::gsl_multifit_nlinear_solver_svd,
            }
        }
    }
}

/// Finite difference type used to approximate the Jacobian.
///
/// Specifies whether to use forward or centered differences when
/// approximating the Jacobian.  This is only used when an analytic
/// Jacobian is not provided to the solver.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FdType {
    /// Specify a forward finite difference to approximate the
    /// Jacobian matrix. The Jacobian matrix will be calculated as
    ///
    /// $$J_{ij} = \frac{1}{Δⱼ\} \bigl( fᵢ(x + Δⱼ eⱼ) - fᵢ(x) \bigr)$$
    ///
    /// where $Δⱼ = h |xⱼ|$ and $eⱼ$ is the standard $j$-th Cartesian
    /// unit basis vector so that $x + Δⱼ eⱼ$ represents a small
    /// (forward) perturbation of the $j$-th parameter by an amount
    /// $Δⱼ$.  The perturbation $Δⱼ$ is proportional to the current
    /// value $|xⱼ|$ which helps to calculate an accurate Jacobian
    /// when the various parameters have different scale sizes.  The
    /// value of h is specified by the h_df parameter.  The accuracy
    /// of this method is $O(h)$, and evaluating this matrix requires
    /// an additional $p$ function evaluations.
    FwDiff,
    /// Specify a centered finite difference to approximate the
    /// Jacobian matrix. The Jacobian matrix will be calculated as
    ///
    /// $$J_{ij} = \frac{1}{Δⱼ} \left(
    /// fᵢ(x + \tfrac{1}{2} Δⱼeⱼ) - fᵢ(x - \tfrac{1}{2} Δⱼeⱼ) \right)$$
    ///
    /// See `FwDiff` for a description of $Δⱼ$.  The accuracy of this
    /// method is $O(h^2)$, but evaluating this matrix requires an
    /// additional $2p$ function evaluations.
    CtrDiff,
}

impl FdType {
    fn to_c(self) -> sys::gsl_multifit_nlinear_fdtype {
        match self {
            Self::FwDiff => sys::gsl_multifit_nlinear_fdtype_GSL_MULTIFIT_NLINEAR_FWDIFF,
            Self::CtrDiff => sys::gsl_multifit_nlinear_fdtype_GSL_MULTIFIT_NLINEAR_CTRDIFF,
        }
    }
}

/// Tunable Parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Parameters {
    /// Trust region subproblem method.
    pub trs: TRS,
    /// Scaling method.
    pub scale: Scale,
    /// Solver method.
    pub solver: Solver,
    /// Finite difference method.
    pub fdtype: FdType,
    /// Factor for increasing trust radius.
    ///
    /// When a step is accepted, the trust region radius will be
    /// increased by this factor.  The default value is 3.
    pub factor_up: f64,
    /// Factor for decreasing trust radius.
    ///
    /// When a step is rejected, the trust region radius will be
    /// decreased by this factor.  The default value is 2.
    pub factor_down: f64,
    /// Max allowed $|a|/|v|$.
    ///
    /// When using geodesic acceleration to solve a nonlinear least
    /// squares problem, an important parameter to monitor is the
    /// ratio of the acceleration term to the velocity term,
    ///
    /// $$\frac{‖a‖}{‖v‖}$$
    ///
    /// If this ratio is small, it means the acceleration correction
    /// is contributing very little to the step.  This could be
    /// because the problem is not “nonlinear” enough to benefit from
    /// the acceleration. If the ratio is large ($> 1$) it means that
    /// the acceleration is larger than the velocity, which shouldn’t
    /// happen since the step represents a truncated series and so the
    /// second order term $a$ should be smaller than the first order
    /// term $v$ to guarantee convergence.  Therefore any steps with a
    /// ratio larger than the parameter `avmax` are rejected.  `avmax`
    /// is set to 0.75 by default.  For problems which experience
    /// difficulty converging, this threshold could be lowered.
    pub avmax: f64,
    /// Step size for finite difference Jacobian.
    ///
    /// This parameter specifies the step size for approximating the
    /// Jacobian matrix with finite differences.  It is set to $√ε$ by
    /// default, where $ε$ is [`f64::EPSILON`].
    pub h_df: f64,
    /// Step size for finite difference $f_{vv}$.
    ///
    /// When using geodesic acceleration, the user must either supply
    /// a function to calculate $f_{vv}(x)$ or the library can
    /// estimate this second directional derivative using a finite
    /// difference method.  When using finite differences, the library
    /// must calculate $f(x + h v)$ where $h$ represents a small step
    /// in the velocity direction.  The parameter `h_fvv` defines this
    /// step size and is set to 0.02 by default.
    pub h_fvv: f64,
}

impl Parameters {
    fn to_c(&self) -> sys::gsl_multifit_nlinear_parameters {
        sys::gsl_multifit_nlinear_parameters {
            trs: self.trs.to_c(),
            scale: self.scale.to_c(),
            solver: self.solver.to_c(),
            fdtype: self.fdtype.to_c(),
            factor_up: self.factor_up,
            factor_down: self.factor_down,
            avmax: self.avmax,
            h_df: self.h_df,
            h_fvv: self.h_fvv,
        }
    }
}

impl Parameters {
    /// Return a set of recommended default parameters for use in
    /// solving nonlinear least squares problems.
    pub fn new() -> Parameters {
        Parameters {
            trs: TRS::LM,
            scale: Scale::More,
            solver: Solver::QR,
            fdtype: FdType::FwDiff,
            factor_up: 3.,
            factor_down: 2.,
            avmax: 0.75,
            h_df: SQRT_EPSILON,
            h_fvv: 0.02,
        }
    }
}

impl std::default::Default for Parameters {
    fn default() -> Self {
        Parameters::new()
    }
}

/// The type of algorithm which will be used to solve a nonlinear
/// least squares problem.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Type {
    /// Trust region method.  It is currently the only implemented
    /// nonlinear least squares method.
    Trust,
}

/// Used to identify the types that can be used to specify the
/// functions to pass to [`Workspace::init`].
pub trait Fdf<V: AsMatrix + ?Sized> {
    /// `f(x, fx)` must store the `n` components of the vector $f(x)$
    /// in `fx` for argument `x`, returning an appropriate error code
    /// if the function cannot be computed.
    fn f(&mut self, x: &V, fx: &mut V) -> Result<(), Error>;

    /// `df(x, J)` stores the `n`-by-`p` matrix result
    ///
    /// $$J_{ij} = ∂fᵢ(x) / ∂xⱼ$$
    ///
    /// in `J` for argument `x`, returning an appropriate error code
    /// if the matrix cannot be computed.  If an analytic Jacobian is
    /// unavailable, or too expensive to compute, you can use
    /// [`unimplemented'()`] and let [`Self::has_df`] return `false`.
    /// In this case the Jacobian will be internally computed using
    /// finite difference approximations of the function f.
    fn df(&mut self, x: &V, J: &mut V::MatViewMut<'_>) -> Result<(), Error>;

    /// Return `true` iff the method [`Self::df`] is implemented.
    fn has_df(&self) -> bool;

    /// `fvv(a, v, fvv)` must store the `n` components of the vector
    /// $f_{vv}(x) = ∑_{αβ} v_α v_β \frac{∂}{∂x_α} \frac{∂}{∂x_β}
    /// f(x)$, representing second directional derivatives of the
    /// function to be minimized, into the `fvv`.  The point is
    /// provided in `x` and the velocity vector is provided in `v`,
    /// both of which have `p` components.
    ///
    /// This is needed when geodesic acceleration is enabled.  If
    /// analytic expressions for $f_{vv}(x)$ are unavailable or too
    /// difficult to compute, this function may be
    /// [`unimplemented!()`] and [`Self::has_fvv`] must return
    /// `false`.  In this case $f_{vv}(x)$ will be computed internally
    /// using a finite difference approximation.
    fn fvv(&mut self, x: &V, v: &V, J: &mut V) -> Result<(), Error>;

    /// Return `true` iff the method [`Self::fvv`] is implemented.
    fn has_fvv(&self) -> bool;
}

impl<V, F> Fdf<V> for F
where
    V: AsMatrix + ?Sized,
    F: FnMut(&V, &mut V) -> Result<(), Error>,
{
    fn f(&mut self, x: &V, fx: &mut V) -> Result<(), Error> {
        self(x, fx)
    }

    #[inline]
    fn has_df(&self) -> bool {
        false
    }

    fn df(&mut self, _x: &V, _J: &mut V::MatViewMut<'_>) -> Result<(), Error> {
        unimplemented!()
    }

    #[inline]
    fn has_fvv(&self) -> bool {
        false
    }

    fn fvv(&mut self, _x: &V, _v: &V, _J: &mut V) -> Result<(), Error> {
        unimplemented!()
    }
}

ffi_wrapper!(
    /// Derivative solver.
    Workspace<'a, V: AsMatrix + ?Sized>,
    *mut sys::gsl_multifit_nlinear_workspace,
    gsl_multifit_nlinear_free
    ;fdf_struct: Option<Box<sys::gsl_multifit_nlinear_fdf>> => None;
    ;fdf: Option<Box<dyn Fdf<V> + 'a>> => None;
);

macro_rules! impl_workspace {
    ($fit: ident, $path: literal) => {
        paste! {

            impl<'a, V: AsMatrix + ?Sized> Workspace<'a, V> {
                fn type_to_c(
                    t: Type
                ) -> *const sys::[<gsl_multi $fit _nlinear_type>] {
                    match t {
                        Type::Trust => unsafe {
                            sys::[<gsl_multi $fit _nlinear_trust>]
                        },
                    }
                }

                /// Return a new instance of a derivative solver of
                /// type `t` for `n` observations and `p` parameters.
                /// The `params` input specifies a tunable set of
                /// parameters which will affect important details in
                /// each iteration of the trust region subproblem
                /// algorithm.  It is recommended to start with the
                /// suggested default parameters
                /// [`Parameters::default()`] and then tune the
                /// parameters once the code is working correctly.
                /// See Tunable [`Parameters`] for descriptions of the
                /// various parameters.  For example, the following
                /// code creates an instance of a Levenberg-Marquardt
                /// solver for `100` data points and `3` parameters,
                /// using suggested defaults:
                ///
                /// ```
                /// use rgsl::VecF64;
                #[doc = "use " $path "::{Parameters, Type, Workspace};"]
                /// let params = Parameters::default();
                /// let w = Workspace::<VecF64>::new(Type::Trust, &params, 100, 3);
                /// ```
                ///
                /// The number of observations `n` must be greater
                /// than or equal to parameters `p`.
                ///
                /// # Panic
                /// Panic if there is insufficient memory to create
                /// the workspace.
                #[doc(alias = gsl_multi $fit _nlinear_alloc)]
                pub fn new(t: Type, params: &Parameters, n: usize, p: usize) -> Self {
                    let w = unsafe { sys::[<gsl_multi $fit _nlinear_alloc>](
                        Self::type_to_c(t), &params.to_c(), n, p)
                    };
                    if w.is_null() {
                        panic!("rgsl::multi{}::Workspace::new: out of memory",
                            stringify!($fit));
                    }
                    Self::wrap(w)
                }

                pub fn trust(params: &Parameters, n: usize, p: usize) -> Self {
                    Self::new(Type::Trust, params, n, p)
                }

                /// Return the number of observations the workspace is for.
                pub fn n(&self) -> usize {
                    unsafe {
                        let w = self.inner.as_ref_unchecked();
                        (*w.f).size
                    }
                }

                /// Return the number of parameters the workspace is for.
                pub fn p(&self) -> usize {
                    unsafe {
                        let w = self.inner.as_ref_unchecked();
                        (*w.x).size
                    }
                }

            }
        }
    };
}

impl_workspace!(fit, "rgsl::multifit");

impl<'a, V: AsMatrix + ?Sized> Workspace<'a, V> {
    /// Initialize, or reinitialize, the workspace to use the function
    /// `f` and the initial guess `x`.
    ///
    /// The call `f(x, fx)` should store the `n` components of the
    /// vector $f(x)$ in `fx` for argument `x`, returning an
    /// appropriate error code if the function cannot be computed.
    ///
    /// You may also want to supply a [function computing the
    /// derivative of $f$](Self::df).
    #[doc(alias = "gsl_multifit_nlinear_init")]
    pub fn init<F: Fdf<V> + 'a>(&mut self, x: &V, fdf: F) -> Result<(), Error> {
        unsafe extern "C" fn f_trampoline<V: AsMatrix + ?Sized, F: Fdf<V>>(
            x: *const sys::gsl_vector,
            params: *mut c_void,
            fx: *mut sys::gsl_vector,
        ) -> i32 {
            let ret = std::panic::catch_unwind(|| unsafe {
                let fdf = &mut *params.cast::<F>();
                let vx = V::view_from_ptr(x);
                let mut vfx = V::view_from_mut_ptr(fx);
                fdf.f(&*vx, &mut *vfx)
            });
            Error::to_c(ret.map_err(|_| Error::Failure).flatten())
        }

        unsafe extern "C" fn df_trampoline<V: AsMatrix + ?Sized, F: Fdf<V>>(
            x: *const sys::gsl_vector,
            params: *mut c_void,
            J: *mut sys::gsl_matrix,
        ) -> i32 {
            let ret = std::panic::catch_unwind(|| unsafe {
                let fdf = &mut *params.cast::<F>();
                let vx = V::view_from_ptr(x);
                let mut vJ = V::mat_view_from_mut_ptr(J);
                fdf.df(&*vx, &mut vJ)
            });
            Error::to_c(ret.map_err(|_| Error::Failure).flatten())
        }

        unsafe extern "C" fn fvv_trampoline<V: AsMatrix + ?Sized, F: Fdf<V>>(
            x: *const sys::gsl_vector,
            v: *const sys::gsl_vector,
            params: *mut c_void,
            fvv: *mut sys::gsl_vector,
        ) -> i32 {
            let ret = std::panic::catch_unwind(|| unsafe {
                let fdf = &mut *params.cast::<F>();
                let vx = V::view_from_ptr(x);
                let vv = V::view_from_ptr(v);
                let mut vfvv = V::view_from_mut_ptr(fvv);
                fdf.fvv(&*vx, &*vv, &mut *vfvv)
            });
            Error::to_c(ret.map_err(|_| Error::Failure).flatten())
        }

        let mut fdf: Box<F> = Box::new(fdf);
        let mut fdf_struct = Box::new(sys::gsl_multifit_nlinear_fdf {
            f: Some(f_trampoline::<V, F>),
            df: fdf.has_df().then_some(df_trampoline::<V, F>),
            fvv: fdf.has_fvv().then_some(fvv_trampoline::<V, F>),
            n: self.n(),
            p: self.p(),
            params: &mut *fdf as *mut F as *mut _,
            nevalf: 0,
            nevaldf: 0,
            nevalfvv: 0,
        });
        self.fdf = Some(fdf);

        let x = V::as_gsl_vector(x);
        let ret = unsafe {
            sys::gsl_multifit_nlinear_init(
                &*x,              // copied into the workspace
                &mut *fdf_struct, // Heap pointer (stable address)
                self.unwrap_unique(),
            )
        };
        self.fdf_struct = Some(fdf_struct);
        Error::handle(ret, ())
    }

    pub fn name(&self) -> TRS {
        let n = unsafe { sys::gsl_multifit_nlinear_trs_name(self.unwrap_shared()) };
        map_name!(
            rgsl::multifit::Workspace,
            [
                (c"levenberg-marquardt", TRS::LM),
                (c"levenberg-marquardt+accel", TRS::LMaccel),
                (c"dogleg", TRS::Dogleg),
                (c"double-dogleg", TRS::DDogleg),
                (c"2D-subspace", TRS::Subspace2D),
            ],
            n,
            TRS
        )
    }
}

pub mod large {
    use super::*;

    /// Method to solve the Trust Region Subproblem.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub enum TRS {
        /// [Levenberg-Marquardt](crate::multifit#levenberg-marquardt)
        /// algorithm.
        LM,
        /// [Levenberg-Marquardt with Geodesic
        /// Acceleration](crate::multifit#levenberg-marquardt-with-geodesic-acceleration).
        LMaccel,
        /// [Dogleg](crate::multifit#dogleg) algorithm.
        Dogleg,
        /// [Double Dogleg](crate::multifit#double-dogleg) algorithm.
        DDogleg,
        /// [Two Dimensional
        /// Subspace](crate::multifit#two-dimensional-subspace) algorithm.
        Subspace2D,
        /// [Steihaug-Toint Conjugate
        /// Gradient](crate::multifit#steihaug-toint-conjugate-gradient)
        /// algorithm. This method is available only for large systems.
        CgST,
    }

    impl TRS {
        fn to_c(self) -> *const sys::gsl_multilarge_nlinear_trs {
            unsafe {
                match self {
                    Self::LM => sys::gsl_multilarge_nlinear_trs_lm,
                    Self::LMaccel => sys::gsl_multilarge_nlinear_trs_lmaccel,
                    Self::Dogleg => sys::gsl_multilarge_nlinear_trs_dogleg,
                    Self::DDogleg => sys::gsl_multilarge_nlinear_trs_ddogleg,
                    Self::Subspace2D => sys::gsl_multilarge_nlinear_trs_subspace2D,
                    Self::CgST => sys::gsl_multilarge_nlinear_trs_cgst,
                }
            }
        }
    }

    /// Make possible to use `TRS` instead of `large::TRS` so switching to
    /// the `large` module requires minimal code change.
    impl std::convert::From<super::TRS> for TRS {
        fn from(trs: super::TRS) -> Self {
            match trs {
                super::TRS::LM => TRS::LM,
                super::TRS::LMaccel => TRS::LMaccel,
                super::TRS::Dogleg => TRS::Dogleg,
                super::TRS::DDogleg => TRS::DDogleg,
                super::TRS::Subspace2D => TRS::Subspace2D,
            }
        }
    }

    pub use super::Scale;

    /// Solver for the trust region subproblem.
    ///
    /// Solving the trust region subproblem on each iteration almost
    /// always requires the solution of the following linear least
    /// squares system
    ///
    /// $$\begin{pmatrix} J \cr √μ D \end{pmatrix} δ
    /// = - \begin{pmatrix} f \cr 0 \end{pmatrix}$$
    ///
    /// How the system is solved and can be selected the choices of this
    /// enum.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub enum Solver {
        /// [Cholesky](super::Solver::Cholesky) method.
        Cholesky,
        /// [Modified Cholesky](super::Solver::MCholesky) method.
        MCholesky,
    }

    impl Solver {
        fn to_c(self) -> *const sys::gsl_multilarge_nlinear_solver {
            unsafe {
                match self {
                    Self::Cholesky => sys::gsl_multilarge_nlinear_solver_cholesky,
                    Self::MCholesky => sys::gsl_multilarge_nlinear_solver_mcholesky,
                }
            }
        }
    }

    pub use super::FdType;

    /// Tunable Parameters.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct Parameters {
        /// Trust region subproblem method.
        pub trs: TRS,
        /// Scaling method.
        pub scale: Scale,
        /// Solver method.
        pub solver: Solver,
        /// Finite difference method.
        pub fdtype: FdType,
        /// Factor for increasing trust radius.
        ///
        /// When a step is accepted, the trust region radius will be
        /// increased by this factor.  The default value is 3.
        pub factor_up: f64,
        /// Factor for decreasing trust radius.
        ///
        /// When a step is rejected, the trust region radius will be
        /// decreased by this factor.  The default value is 2.
        pub factor_down: f64,
        /// Max allowed $|a|/|v|$.
        ///
        /// See [`Parameters::avmax`](super::Parameters::avmax) for
        /// more details.
        pub avmax: f64,
        /// Step size for finite difference Jacobian.
        ///
        /// See [`Parameters::h_df`](super::Parameters::h_df) for more
        /// details.
        pub h_df: f64,
        /// Step size for finite difference $f_{vv}$.
        ///
        /// See [`Parameters::h_fvv`](super::Parameters::h_fvv) for
        /// more details.
        pub h_fvv: f64,
        /// Maximum iterations for trs method.
        pub max_iter: usize,
        /// Tolerance for solving trs.
        pub tol: f64,
    }

    impl Parameters {
        fn to_c(&self) -> sys::gsl_multilarge_nlinear_parameters {
            sys::gsl_multilarge_nlinear_parameters {
                trs: self.trs.to_c(),
                scale: self.scale.to_c_large(),
                solver: self.solver.to_c(),
                fdtype: self.fdtype.to_c(),
                factor_up: self.factor_up,
                factor_down: self.factor_down,
                avmax: self.avmax,
                h_df: self.h_df,
                h_fvv: self.h_fvv,
                max_iter: self.max_iter,
                tol: self.tol,
            }
        }

        /// Return the parameters set to their default values.
        pub fn new() -> Self {
            Self {
                trs: TRS::LM,
                scale: Scale::More,
                solver: Solver::Cholesky,
                fdtype: FdType::FwDiff,
                factor_up: 3.,
                factor_down: 2.,
                avmax: 0.75,
                h_df: SQRT_EPSILON,
                h_fvv: 0.01,
                max_iter: 0,
                tol: 1e-6,
            }
        }
    }

    impl std::default::Default for Parameters {
        fn default() -> Self {
            Parameters::new()
        }
    }

    pub use super::Type;

    /// Used to identify the types that can be used to specify the
    /// functions to pass to [`Workspace::init`].
    pub trait Fdf<V: AsMatrix + ?Sized> {
        /// `f(x, fx)` must store the `n` components of the vector $f(x)$
        /// in `fx` for argument `x`, returning an appropriate error code
        /// if the function cannot be computed.
        fn f(&mut self, x: &V, fx: &mut V) -> Result<(), Error>;

        /// If `transJ` is `false`, then this function should compute
        /// the matrix-vector product $J u$, where
        ///
        /// $$J_{ij} = ∂fᵢ(x) / ∂xⱼ$$
        ///
        /// and store the result in `v`.  If `transJ` is `true`, then
        /// this function should compute the matrix-vector product
        /// $J^T u$ and store the result in `v`.  Additionally, the
        /// normal equations matrix $J^T J$ should be stored in the
        /// lower half of `JTJ`.  The matrix `JTJ` may not be
        /// provided, for example by iterative methods which do not
        /// require this matrix, so the user should check for this
        /// prior to constructing the matrix.
        fn df(
            &mut self,
            transJ: bool,
            x: &V,
            u: &V,
            v: &mut V,
            JTJ: Option<&mut V::MatViewMut<'_>>,
        ) -> Result<(), Error>;

        /// `fvv(a, v, fvv)` must store the `n` components of the vector
        /// $f_{vv}(x) = ∑_{αβ} v_α v_β \frac{∂}{∂x_α} \frac{∂}{∂x_β}
        /// f(x)$, representing second directional derivatives of the
        /// function to be minimized, into the `fvv`.  The point is
        /// provided in `x` and the velocity vector is provided in `v`,
        /// both of which have `p` components.
        ///
        /// This is needed when geodesic acceleration is enabled.  If
        /// analytic expressions for $f_{vv}(x)$ are unavailable or too
        /// difficult to compute, this function may be
        /// [`unimplemented!()`] and [`Self::has_fvv`] must return
        /// `false`.  In this case $f_{vv}(x)$ will be computed internally
        /// using a finite difference approximation.
        fn fvv(&mut self, x: &V, v: &V, J: &mut V) -> Result<(), Error>;

        /// Return `true` iff the method [`Self::fvv`] is implemented.
        fn has_fvv(&self) -> bool;
    }

    ffi_wrapper!(
        /// Derivative solver for large problems.
        Workspace<'a, V: AsMatrix + ?Sized>,
        *mut sys::gsl_multilarge_nlinear_workspace,
        gsl_multilarge_nlinear_free
        ;fdf_struct: Option<Box<sys::gsl_multilarge_nlinear_fdf>> => None;
        ;fdf: Option<Box<dyn Fdf<V> + 'a>> => None;
    );

    impl_workspace!(large, "rgsl::multilarge");

    impl<'a, V: AsMatrix + ?Sized> Workspace<'a, V> {
        /// Initialize, or reinitialize, the workspace to use the function
        /// `f` and the initial guess `x`.
        ///
        /// The call `f(x, fx)` should store the `n` components of the
        /// vector $f(x)$ in `fx` for argument `x`, returning an
        /// appropriate error code if the function cannot be computed.
        ///
        /// You may also want to supply a [function computing the
        /// derivative of $f$](Self::df).
        #[doc(alias = "gsl_multilarge_nlinear_init")]
        pub fn init<F: Fdf<V> + 'a>(&mut self, x: &V, fdf: F) -> Result<(), Error> {
            unsafe extern "C" fn f_trampoline<V: AsMatrix + ?Sized, F: Fdf<V>>(
                x: *const sys::gsl_vector,
                params: *mut c_void,
                fx: *mut sys::gsl_vector,
            ) -> i32 {
                let ret = std::panic::catch_unwind(|| unsafe {
                    let fdf = &mut *params.cast::<F>();
                    let vx = V::view_from_ptr(x);
                    let mut vfx = V::view_from_mut_ptr(fx);
                    fdf.f(&*vx, &mut *vfx)
                });
                Error::to_c(ret.map_err(|_| Error::Failure).flatten())
            }

            unsafe extern "C" fn df_trampoline<V: AsMatrix + ?Sized, F: Fdf<V>>(
                transJ: sys::CBLAS_TRANSPOSE,
                x: *const sys::gsl_vector,
                u: *const sys::gsl_vector,
                params: *mut c_void,
                v: *mut sys::gsl_vector,
                JTJ: *mut sys::gsl_matrix,
            ) -> i32 {
                let transJ = transJ == sys::CBLAS_TRANSPOSE_CblasTrans;
                let ret = std::panic::catch_unwind(|| unsafe {
                    let fdf = &mut *params.cast::<F>();
                    let vx = V::view_from_ptr(x);
                    let vu = V::view_from_ptr(u);
                    let mut vv = V::view_from_mut_ptr(v);
                    let mut vJTJ = if JTJ.is_null() {
                        None
                    } else {
                        Some(V::mat_view_from_mut_ptr(JTJ))
                    };
                    fdf.df(transJ, &*vx, &*vu, &mut *vv, vJTJ.as_mut())
                });
                Error::to_c(ret.map_err(|_| Error::Failure).flatten())
            }

            unsafe extern "C" fn fvv_trampoline<V: AsMatrix + ?Sized, F: Fdf<V>>(
                x: *const sys::gsl_vector,
                v: *const sys::gsl_vector,
                params: *mut c_void,
                fvv: *mut sys::gsl_vector,
            ) -> i32 {
                let ret = std::panic::catch_unwind(|| unsafe {
                    let fdf = &mut *params.cast::<F>();
                    let vx = V::view_from_ptr(x);
                    let vv = V::view_from_ptr(v);
                    let mut vfvv = V::view_from_mut_ptr(fvv);
                    fdf.fvv(&*vx, &*vv, &mut *vfvv)
                });
                Error::to_c(ret.map_err(|_| Error::Failure).flatten())
            }

            let mut fdf: Box<F> = Box::new(fdf);
            let mut fdf_struct = Box::new(sys::gsl_multilarge_nlinear_fdf {
                f: Some(f_trampoline::<V, F>),
                df: Some(df_trampoline::<V, F>),
                fvv: fdf.has_fvv().then_some(fvv_trampoline::<V, F>),
                n: self.n(),
                p: self.p(),
                params: &mut *fdf as *mut F as *mut _,
                nevalf: 0,
                nevaldfu: 0,
                nevaldf2: 0,
                nevalfvv: 0,
            });
            self.fdf = Some(fdf);

        let x = V::as_gsl_vector(x);
        let ret = unsafe {
            sys::gsl_multilarge_nlinear_init(
                &*x,              // copied into the workspace
                &mut *fdf_struct, // Heap pointer (stable address)
                self.unwrap_unique(),
            )
        };
        self.fdf_struct = Some(fdf_struct);
        Error::handle(ret, ())
    }

        pub fn name(&self) -> TRS {
            let n = unsafe { sys::gsl_multilarge_nlinear_trs_name(self.unwrap_shared()) };
            map_name!(
                rgsl::multifit::Workspace,
                [
                    (c"levenberg-marquardt", TRS::LM),
                    (c"levenberg-marquardt+accel", TRS::LMaccel),
                    (c"dogleg", TRS::Dogleg),
                    (c"double-dogleg", TRS::DDogleg),
                    (c"2D-subspace", TRS::Subspace2D),
                    (c"steihaug-toint", TRS::CgST),
                ],
                n,
                TRS
            )
        }
    }

    // FIXME: keep ?
    #[doc(alias = "gsl_multilarge_linear_L_decomp")]
    pub fn linear_L_decomp(L: &mut MatF64, tau: &mut VecF64) -> Result<(), Error> {
        let ret =
            unsafe { sys::gsl_multilarge_linear_L_decomp(L.unwrap_unique(), tau.unwrap_unique()) };
        Error::handle(ret, ())
    }
}

/// Compute the covariance matrix cov = inv (J^T J) by QRP^T decomposition of J
#[doc(alias = "gsl_multifit_covar")]
pub fn covar(J: &MatF64, epsrel: f64, covar: &mut MatF64) -> Result<(), Error>
 {
    let ret = unsafe { sys::gsl_multifit_covar(J.unwrap_shared(), epsrel, covar.unwrap_unique()) };
    Error::handle(ret, ())
}

#[doc(alias = "gsl_multifit_test_delta")]
pub fn test_delta(dx: &VecF64, x: &VecF64, epsabs: f64, epsrel: f64) -> Result<(), Error> {
    let ret = unsafe {
        sys::gsl_multifit_test_delta(dx.unwrap_shared(), x.unwrap_shared(), epsabs, epsrel)
    };
    Error::handle(ret, ())
}

#[doc(alias = "gsl_multifit_gradient")]
pub fn gradient(J: &MatF64, f: &VecF64, g: &mut VecF64) -> Result<(), Error> {
    let ret = unsafe {
        sys::gsl_multifit_gradient(J.unwrap_shared(), f.unwrap_shared(), g.unwrap_unique())
    };
    Error::handle(ret, ())
}

#[doc(alias = "gsl_multifit_linear_lreg")]
pub fn linear_lreg(smin: f64, smax: f64, reg_param: &mut VecF64) -> Result<(), Error> {
    let ret = unsafe { sys::gsl_multifit_linear_lreg(smin, smax, reg_param.unwrap_unique()) };
    Error::handle(ret, ())
}

/// Returns `idx`.
#[doc(alias = "gsl_multifit_linear_lcorner")]
pub fn linear_lcorner(rho: &VecF64, eta: &VecF64) -> Result<usize, Error> {
    let mut idx = 0;
    let ret = unsafe {
        sys::gsl_multifit_linear_lcorner(rho.unwrap_shared(), eta.unwrap_shared(), &mut idx)
    };
    Error::handle(ret, idx)
}

/// Returns `(Value, idx)`.
#[doc(alias = "gsl_multifit_linear_lcorner2")]
pub fn linear_lcorner2(rho: &VecF64, eta: &VecF64) -> Result<usize, Error> {
    let mut idx = 0;
    let ret = unsafe {
        sys::gsl_multifit_linear_lcorner2(rho.unwrap_shared(), eta.unwrap_shared(), &mut idx)
    };
    Error::handle(ret, idx)
}

#[doc(alias = "gsl_multifit_linear_Lk")]
pub fn linear_Lk(p: usize, k: usize, L: &mut MatF64) -> Result<(), Error> {
    let ret = unsafe { sys::gsl_multifit_linear_Lk(p, k, L.unwrap_unique()) };
    Error::handle(ret, ())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_name() {
        type W<'a> = Workspace<'a, VecF64>;

        let params = Parameters::default();
        let w = W::trust(&params, 10, 2);
        assert_eq!(w.name(), TRS::LM);

        let p = Parameters {
            trs: TRS::LMaccel,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).name(), TRS::LMaccel);

        let p = Parameters {
            trs: TRS::Dogleg,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).name(), TRS::Dogleg);

        let p = Parameters {
            trs: TRS::DDogleg,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).name(), TRS::DDogleg);

        let p = Parameters {
            trs: TRS::Subspace2D,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).name(), TRS::Subspace2D);
    }

    #[test]
    fn test_large_name() {
        use large::{Parameters, TRS};
        type W<'a> = large::Workspace<'a, VecF64>;

        let params = Parameters::default();
        let w = W::trust(&params, 10, 2);
        assert_eq!(w.name(), TRS::LM);

        let p = Parameters {
            trs: TRS::LMaccel,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).name(), TRS::LMaccel);

        let p = Parameters {
            trs: TRS::Dogleg,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).name(), TRS::Dogleg);

        let p = Parameters {
            trs: TRS::DDogleg,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).name(), TRS::DDogleg);

        let p = Parameters {
            trs: TRS::Subspace2D,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).name(), TRS::Subspace2D);

            let p = Parameters {
            trs: TRS::CgST,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).name(), TRS::CgST);
}
}
