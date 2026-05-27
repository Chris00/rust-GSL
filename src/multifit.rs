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


/// Compute the covariance matrix cov = inv (J^T J) by QRP^T decomposition of J
#[doc(alias = "gsl_multifit_covar")]
pub fn covar(J: &MatrixF64, epsrel: f64, covar: &mut MatrixF64) -> Result<(), Error> {
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
