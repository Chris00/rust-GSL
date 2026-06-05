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
to moderate sized problems ([`Workspace`] interface), this is
accomplished by factoring the full Jacobian matrix, which is provided
by the user, with the Cholesky, QR, or SVD decompositions.  For large
systems ([`large::Workspace`] interface), the user has two
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
residuals and Jacobian, or use the [`Workspace::winit`]
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

## Troubleshooting

When developing a code to solve a nonlinear least squares problem,
here are a few considerations to keep in mind.

1. The most common difficulty is the accurate implementation of the
   Jacobian matrix.  If the analytic Jacobian is not properly provided
   to the solver, this can hinder and many times prevent convergence
   of the method.  When developing a new nonlinear least squares code,
   it often helps to compare the program output with the internally
   computed finite difference Jacobian and the user supplied analytic
   Jacobian.  If there is a large difference in coefficients, it is
   likely the analytic Jacobian is incorrectly implemented.
2. If your code is having difficulty converging, the next thing to
   check is the starting point provided to the solver.  The methods of
   this chapter are local methods, meaning if you provide a starting
   point far away from the true minimum, the method may converge to a
   local minimum or not converge at all.  Sometimes it is possible to
   solve a linearized approximation to the nonlinear problem, and use
   the linear solution as the starting point to the nonlinear problem.

3. If the various parameters of the coefficient vector $x$ vary widely
   in magnitude, then the problem is said to be badly scaled.  The
   methods of this chapter do attempt to automatically rescale the
   elements of $x$ to have roughly the same order of magnitude, but in
   extreme cases this could still cause problems for convergence.  In
   these cases it is recommended for the user to scale their parameter
   vector $x$ so that each parameter spans roughly the same range, say
   $[-1,1]$.  The solution vector can be backscaled to recover the
   original units of the problem.

## Examples

The following example programs demonstrate the nonlinear least squares
fitting capabilities.

### Exponential Fitting Example

The following example program fits a weighted exponential model with
background to experimental data, $Y = A \exp(-λ t) + b$.  The first
part of the program sets up the functions `expb_f()` and `expb_df()`
to calculate the model and its Jacobian.  The appropriate fitting
function is given by,

$$fᵢ(A, λ, b) := (A \exp(-λtᵢ) + b) - yᵢ$$

where we have chosen $tᵢ := i T / (N - 1)$, where $N$ is the number of
data points fitted, so that $tᵢ ∈ [0, T]$.  The Jacobian matrix $J$ is
the derivative of these functions with respect to the three parameters
$(A, λ, b)$.  It is given by,

$$J_{ij} = \frac{∂fᵢ}{∂xⱼ}$$

where $x₀ = A$, $x₁ = λ$ and $x₂ = b$.  The $i$-th row of the Jacobian
is therefore

$$J_{i\cdot} = \begin{pmatrix} \exp(-λtᵢ) & -tᵢ A \exp(-λtᵢ) & 1 \end{pmatrix}$$

The main part of the program sets up a Levenberg-Marquardt solver and
some simulated random data.  The data uses the known parameters $(5.0,
1.5, 1.0)$ combined with Gaussian noise (standard deviation = 0.1)
with a maximum time $T = 3$ and $N = 100$ timesteps.  The initial
guess for the parameters is chosen as $(1.0, 1.0, 0.0)$.  The
iteration terminates when the relative change in $x$ is smaller than
$10^{-8}$, or when the magnitude of the gradient falls below
$10^{-8}$. Here are the results of running the program:

```text
iter:  0: a = 1.0000, λ = 1.0000, b = 0.0000, cond(J) =      inf, |f(x)| = 88.4448
iter:  1: a = 4.5109, λ = 2.5258, b = 1.0704, cond(J) =  26.2686, |f(x)| = 24.0646
iter:  2: a = 4.8565, λ = 1.7442, b = 1.1669, cond(J) =  23.7470, |f(x)| = 11.9797
iter:  3: a = 4.9356, λ = 1.5713, b = 1.0767, cond(J) =  17.5849, |f(x)| = 10.7355
iter:  4: a = 4.8678, λ = 1.4838, b = 1.0252, cond(J) =  16.3428, |f(x)| = 10.5000
iter:  5: a = 4.8118, λ = 1.4481, b = 1.0076, cond(J) =  15.7925, |f(x)| = 10.4786
iter:  6: a = 4.7983, λ = 1.4404, b = 1.0041, cond(J) =  15.5840, |f(x)| = 10.4778
iter:  7: a = 4.7967, λ = 1.4395, b = 1.0037, cond(J) =  15.5396, |f(x)| = 10.4778
iter:  8: a = 4.7965, λ = 1.4394, b = 1.0037, cond(J) =  15.5344, |f(x)| = 10.4778
iter:  9: a = 4.7965, λ = 1.4394, b = 1.0037, cond(J) =  15.5339, |f(x)| = 10.4778
iter: 10: a = 4.7965, λ = 1.4394, b = 1.0037, cond(J) =  15.5339, |f(x)| = 10.4778
iter: 11: a = 4.7965, λ = 1.4394, b = 1.0037, cond(J) =  15.5339, |f(x)| = 10.4778
Summary from method Trust/LM
Number of iterations: 11
Function evaluations: 16
Jacobian evaluations: 12
Reason for stopping: Gradient
initial |f(x)| = 88.444756
final   |f(x)| = 10.477801
χ²/dof = 1.1318
a = 4.79653 ± 0.18704
λ = 1.43937 ± 0.07390
b = 1.00368 ± 0.03473
```

The approximate values of the parameters are found correctly, and the
chi-squared value indicates a good fit (the chi-squared per degree of
freedom is approximately 1).  In this case the errors on the
parameters can be estimated from the square roots of the diagonal
elements of the covariance matrix.  If the chi-squared value shows a
poor fit (i.e. $χ²/(n-p) ≫ 1$ then the error estimates obtained
from the covariance matrix will be too small. In the example program
the error estimates are multiplied by $√{χ²/(n-p)}$ in this
case, a common way of increasing the errors for a poor fit.  Note that
a poor fit will result from the use of an inappropriate model, and the
scaled error estimates may then be outside the range of validity for
Gaussian errors.

Additionally, we see that the condition number of $J(x)$ stays
reasonably small throughout the iteration.  This indicates we could
safely switch to the Cholesky solver for speed improvement, although
this particular system is too small to really benefit.

![Exponential fitted curve with data](https://www.gnu.org/software/gsl/doc/html/_images/fit-exp.png)

```
use rgsl::{
    Rng, RngType,
    MatF64, VecF64, Error,
    blas::d::{dot, nrm2},
    multifit::{self, Parameters},
};

const TMAX: f64 = 3.; // time variable in [0,TMAX]

struct Data {
    t: VecF64,
    y: VecF64,
    weights: VecF64,
}

impl Data {
    /// Return the sampling data with `n` points.
    fn new(n: usize) -> Self {
        let mut r = Rng::new(RngType::default());
        let mut t = VecF64::zeros(n);
        let mut y = VecF64::zeros(n);
        let mut weights = VecF64::zeros(n);
        for i in 0..n {
            t[i] = i as f64 * TMAX / (n - 1) as f64;
            let yi = 1. + 5. * (-1.5 * t[i]).exp();
            let σi = 0.1 * yi;
            let dy = r.gaussian(σi);
            y[i] = yi + dy;
            weights[i] = 1. / σi.powi(2);
            println!("data: {} {} {}", t[i], y[i], σi);
        }
        Self { t, y, weights }
    }

    fn expb_f(&self, x: &VecF64, fx: &mut VecF64) -> Result<(), Error> {
        let (a, λ, b) = (x[0], x[1], x[2]);
        for i in 0..fx.len() {
            // Model Yᵢ = a exp(-λtᵢ) + b
            let yi = a * (-λ * self.t[i]).exp() + b;
            fx[i] = yi - self.y[i];
        }
        Ok(())
    }

    fn expb_df(&self, x: &VecF64, jac: &mut MatF64) -> Result<(), Error> {
        let (a, λ) = (x[0], x[1]);
        for i in 0..jac.nrows() {
            // Jacobian matrix Jᵢⱼ = ∂fᵢ/∂xⱼ,
            // where fᵢ = (Yᵢ - yᵢ)/σ[i],
            //       Yᵢ = a exp(-λ tᵢ) + b
            // and the xⱼ are the parameters (a,λ,b)
            let e = (-λ * self.t[i]).exp();
            jac[(i, 0)] = e;
            jac[(i, 1)] = -self.t[i] * a * e;
            jac[(i, 2)] = 1.;
        }
        Ok(())
    }
}

fn main() -> Result<(), Error> {
    let params = Parameters::default();
    let d = Data::new(100);
    let x: VecF64 = [1., 1., 0.].into();
    let mut w = multifit::Workspace::trust(&params, d.t.len(), 3)
        .winit(&x, &d.weights,
        (|x: &_, fx: &mut _| d.expb_f(x, fx),
         |x: &_, jac: &mut _| d.expb_df(x, jac)))?;

    // Compute initial cost function.
    let chisq0 = {
        let f = w.residual();
        dot(&f, &f)?
    };
    let info = w.driver(100)
        .xtol(1e-8).gtol(1e-8).ftol(0.)
        .cb(|w, iter| {
            let x = w.position();
            let fx = w.residual();
            let rcond = w.rcond().unwrap();
            println!("iter: {:>2}: a = {:.4}, λ = {:.4}, b = {:.4}, \
                cond(J) = {:8.4}, |f(x)| = {:.4}",
                iter, x[0], x[1], x[2], 1. / rcond, nrm2(&fx));
        })
        .run()?;
    let jac = w.jac();
    let mut covar = MatF64::zeros(x.len(), x.len());
    multifit::covar(&jac, 0., &mut covar)?;
    // Compute final cost.
    let chisq = {
        let f = w.residual();
        dot(&f, &f)?
    };
    println!("Summary from method {:?}/{:?}", w.name(), w.trs_name());
    println!("Number of iterations: {}", w.niter());
    println!("Function evaluations: {}", w.nevalf());
    println!("Jacobian evaluations: {}", w.nevaldf());
    println!("Reason for stopping: {:?}", info);
    println!("initial |f(x)| = {:.6}", chisq0.sqrt());
    println!("final   |f(x)| = {:.6}", chisq.sqrt());
    let dof = (w.n() - w.p()) as f64;
    let c = (chisq / dof).sqrt().max(1.);
    println!("χ²/dof = {:.4}", chisq / dof);
    let x = w.position();
    for (l, i) in [("a", 0), ("λ", 1), ("b", 2)] {
        println!("{l} = {:.5} ± {:.5}", x[i], c * covar[(i,i)].sqrt());
    }
    # // assert!(false); // Uncomment to see output
    Ok(())
}
```

### Geodesic Acceleration Example 1

The following example program minimizes a modified Rosenbrock
function, which is characterized by a narrow canyon with steep
walls. The starting point is selected high on the canyon wall, so the
solver must first find the canyon bottom and then navigate to the
minimum. The problem is solved both with and without using geodesic
acceleration for comparison. The cost function is given by

$$\begin{aligned}
Φ(x) &= \frac{1}{2} (f₁² + f₂²) \cr
f₁ &= 100 (x₂ - x₁²) \cr
f₂ &= 1 - x_1
\end{aligned}$$

The Jacobian matrix is

$$J = \begin{pmatrix}
∂f₁/∂x₁ & ∂f₁/∂x₂ \cr
∂f₂/∂x₁ & ∂f₂/∂x₂
\end{pmatrix}
= \begin{pmatrix} -200 x_1 & 100 \cr -1 & 0 \end{pmatrix}.$$

In order to use geodesic acceleration, the user must provide the
second directional derivative of each residual in the velocity
direction, $Dᵥ² fᵢ = ∑_{αβ} v_α v_β ∂_α ∂_β fᵢ$.  The velocity vector
$v$ is provided by the solver.  For this example, these derivatives
are

$$f_{vv} = Dᵥ² \begin{pmatrix} f₁ \cr f₂ \end{pmatrix}
= \begin{pmatrix} -200 v₁² \cr 0 \end{pmatrix}.$$

The solution of this minimization problem is

$$\begin{aligned}
x^{\*} &= \begin{pmatrix} 1 \cr 1 \end{pmatrix}\cr
Φ(x^{\*}) &= 0
\end{aligned}$$

The program output is shown below:

```text
=== Solving system without acceleration ===
NITER         = 53
NFEV          = 56
NJEV          = 54
NAEV          = 0
initial cost  = 2.250225000000e4
final cost    = 6.674986031430e-18
final x       = [0.9999999974164856, 0.9999999948327616]
final cond(J) = 6.000096055094e2

=== Solving system with acceleration ===
NITER         = 15
NFEV          = 17
NJEV          = 16
NAEV          = 16
initial cost  = 2.250225000000e4
final cost    = 7.518932873279e-19
final x       = [0.999999999132885, 0.999999998265748]
final cond(J) = 6.000097233278e2
```

![Paths taken by solver for Rosenbrock
 function](https://www.gnu.org/software/gsl/doc/html/_images/nlfit2.png)

We can see that enabling geodesic acceleration requires less than a
third of the number of Jacobian evaluations in order to locate the
minimum. The path taken by both methods is shown in the above figure.
The contours show the cost function $Φ(x₁,x₂)$.  We see that both
methods quickly find the canyon bottom, but the geodesic acceleration
method navigates along the bottom to the solution with significantly
fewer iterations.

The program is given below.

```
use rgsl::{
    MatF64, VecF64, Error,
    blas::d::dot,
    multifit::{Fdf, Parameters, TRS, Workspace},
};

fn f(x: &VecF64, fx: &mut VecF64) -> Result<(), Error> {
    let (x1, x2) = (x[0], x[1]);
    fx[0] = 100. * (x2 - x1.powi(2));
    fx[1] = 1. - x1;
    Ok(())
}

fn df(x: &VecF64, jac: &mut MatF64) -> Result<(), Error> {
    let x1 = x[0];
    jac[(0, 0)] = -200. * x1;
    jac[(0, 1)] = 100.;
    jac[(1, 0)] = -1.;
    jac[(1, 1)] = 0.;
    Ok(())
}

fn fvv(x: &VecF64, v: &VecF64, fvv: &mut VecF64) -> Result<(), Error> {
    fvv[0] = -200. * v[0].powi(2);
    fvv[1] = 0.;
    Ok(())
}

fn callback(w: &Workspace<VecF64>, iter: usize) {
    eprintln!("{:?}", w.position());
}

fn solve_system(
    x0: &VecF64,
    n: usize,
    fdf: impl Fdf<VecF64>,
    params: &Parameters,
) -> Result<(), Error> {
    let max_iter = 200;
    let mut work = Workspace::trust(params, n, x0.len())
        .init(x0, fdf)?;
    let χsq0 = {
        let f = work.residual();
        dot(&f, &f)?
    };
    // Iterate until convergence.
    work.driver(max_iter)
        .xtol(1e-8).gtol(1e-8).ftol(1e-8)
        .cb(callback)
        .run()?;
    let χsq = {
        let f = work.residual();
        dot(&f, &f)?
    };
    // Store cond(J(x)).
    let rcond = work.rcond().unwrap();
    println!("NITER         = {}", work.niter());
    println!("NFEV          = {}", work.nevalf());
    println!("NJEV          = {}", work.nevaldf());
    println!("NAEV          = {}", work.nevalfvv());
    println!("initial cost  = {:.12e}", χsq0);
    println!("final cost    = {:.12e}", χsq);
    println!("final x       = {:?}", work.position());
    println!("final cond(J) = {:.12e}", 1. / rcond);
    Ok(())
}

fn main() -> Result<(), Error> {
    let n = 2;
    let mut params = Parameters::default();
    let x: VecF64 = [-0.5, 1.75].into();
    println!("=== Solving system without acceleration ===");
    params.trs = TRS::LM;
    solve_system(&x, n, (f, df, fvv), &params)?;

    println!("\n=== Solving system with acceleration ===");
    params.trs = TRS::LMaccel;
    solve_system(&x, n, (f, df, fvv), &params)?;

    # // assert!(false); // Uncomment to see output
    Ok(())
}
```

### Geodesic Acceleration Example 2

The following example fits a set of data to a Gaussian model using the
Levenberg-Marquardt method with geodesic acceleration. The cost
function is

$$\begin{aligned}
Φ(x) &= \frac{1}{2} ∑_{i=1}^n fᵢ² \cr
fᵢ &= yᵢ - Y(a,b,c; tᵢ)
\end{aligned}$$

where $x = (a, b, c)$ and $yᵢ$ is the measured data point at time
$tᵢ$, and the model is specified by

$$Y(a,b,c; t)
= a \exp\left[ -\frac{1}{2} \left(\frac{t - b}{c} \right)² \right].$$

The parameters $a,b,c$ represent the amplitude, mean, and width of the
Gaussian respectively.  The program below generates the $yᵢ$ data on
$\[0,1\]$ using the values $a = 5$, $b = 0.4$, $c = 0.15$ and adding
random noise.  The $i$-th row of the Jacobian is

$$J_{i,:} = \begin{pmatrix}
\frac{∂fᵢ}{∂a} & \frac{∂fᵢ}{∂b} & \frac{∂fᵢ}{∂c}
\end{pmatrix}
= \begin{pmatrix}
-eᵢ & -\frac{a}{c} zᵢ eᵢ & -\frac{a}{c} zᵢ² eᵢ
\end{pmatrix}$$

where

$$\begin{aligned}
zᵢ &= \frac{tᵢ - b}{c} \cr
eᵢ &= \exp\left( -\frac{1}{2} zᵢ² \right)
\end{aligned}$$

In order to use geodesic acceleration, we need the second directional
derivative of the residuals in the velocity direction, $Dᵥ² fᵢ =
∑_{αβ} v_α v_β ∂_α ∂_β fᵢ$, where $v$ is provided by the solver.  To
compute this, it is helpful to make a table of all second derivatives
of the residuals $fᵢ$ with respect to each combination of model
parameters.  This table is

$$\begin{array}{cccc}
            & \frac{∂}{∂a} & \frac{∂}{∂b} & \frac{∂}{∂c} \cr
\frac{∂}{∂a}& 0 & -\frac{zᵢ}{c} eᵢ & -\frac{zᵢ²}{c} eᵢ\cr
\frac{∂}{∂b}& & \frac{a}{c²} (1 - zᵢ²) eᵢ & \frac{a}{c²} zᵢ (2 - zᵢ²) eᵢ\cr
\frac{∂}{∂c}& & & \frac{a}{c²} zᵢ²(3 - zᵢ²) eᵢ
\end{array}$$

The lower half of the table is omitted since it is symmetric.  Then,
the second directional derivative of $fᵢ$ is

$$Dᵥ² fᵢ
= vₐ² ∂ₐ² fᵢ + 2 vₐ v\_b ∂ₐ ∂\_b fᵢ + 2 vₐ v\_c ∂ₐ ∂\_c fᵢ +
v\_b² ∂\_b² fᵢ + 2 v\_b v\_c ∂\_b ∂\_c fᵢ + v\_c² ∂\_c² fᵢ.$$

The factors of 2 come from the symmetry of the mixed second partial
derivatives.  The iteration is started using the initial guess $a =
1$, $b = 0$, $c = 1$.  The program output is shown below:

```text
            a      b      c    |a|/|v|  cond(J)  |f(x)|
iter  0: 1.0000 0.0000 1.0000   0.0000      inf  35.4785
iter  1: 1.5708 0.5321 0.5219   0.3093  29.0443  31.1042
iter  2: 1.7387 0.4040 0.4568   0.1199   3.5256  28.7217
iter  3: 2.2340 0.3829 0.3053   0.3308   4.5121  23.8074
iter  4: 3.2275 0.3952 0.2243   0.2784   8.6499  15.6003
iter  5: 4.3347 0.3974 0.1752   0.2029  15.1732   7.5908
iter  6: 4.9352 0.3992 0.1536   0.1001  26.6621   4.8402
iter  7: 5.0716 0.3994 0.1498   0.0166  34.6922   4.7103
iter  8: 5.0828 0.3994 0.1495   0.0012  36.5422   4.7095
iter  9: 5.0831 0.3994 0.1495   0.0000  36.6929   4.7095
iter 10: 5.0831 0.3994 0.1495   0.0000  36.6975   4.7095
iter 11: 5.0831 0.3994 0.1495   0.0000  36.6976   4.7095
NITER         = 11
NFEV          = 18
NJEV          = 12
NAEV          = 17
initial cost  = 1.258724737288e3
final cost    = 2.217977560180e1
final x       = [5.083101559156272, 0.39944841095936445, 0.14948983066065347]
final cond(J) = 3.669757713403e1
```

We see the method converges after 11 iterations.  For comparison the
standard Levenberg-Marquardt method requires 26 iterations and so the
Gaussian fitting problem benefits substantially from the geodesic
acceleration correction.  The column marked $|a|/|v|$ above shows the
ratio of the acceleration term to the velocity term as the iteration
progresses.  Larger values of this ratio indicate that the geodesic
acceleration correction term is contributing substantial information
to the solver relative to the standard LM velocity step.

The data and fitted model are shown in the next figure.

![Gaussian model fitted to data](https://www.gnu.org/software/gsl/doc/html/_images/nlfit2b.png)

```
use rgsl::{
    Error, MatF64, VecF64,
    blas::d::{dot, nrm2},
    rng::{Rng, RngType},
    multifit::{Fdf, TRS, Workspace, Parameters},
};

// Model function: a exp(-1/2 [(t - b)/c]²)
fn gaussian(a: f64, b: f64, c: f64, t: f64) -> f64 {
    let z = (t - b) / c;
    a * (-0.5 * z * z).exp()
}

struct Data {
    t: VecF64,
    y: VecF64,
}

impl Data {
    // Generate synthetic data with noise with `n` points.
    fn new(a: f64, b: f64, c: f64, n: usize) -> Self {
        let mut r = Rng::new(RngType::default());
        let mut t = VecF64::zeros(n);
        let mut y = VecF64::zeros(n);
        for i in 0..n {
            t[i] = i as f64 / n as f64;
            let y0 = gaussian(a, b, c, t[i]);
            let dy = r.gaussian(0.1 * y0);
            y[i] = y0 + dy;
        }
        Self { t, y }
    }

    fn f(&self, x: &VecF64, fx: &mut VecF64) -> Result<(), Error> {
        let (a, b, c) = (x[0], x[1], x[2]);
        for i in 0..self.t.len() {
            let y = gaussian(a, b, c, self.t[i]);
            fx[i] = self.y[i] - y;
        }
        Ok(())
    }

    fn df(&self, x: &VecF64, jac: &mut MatF64) -> Result<(), Error> {
        let (a, b, c) = (x[0], x[1], x[2]);
        for i in 0..self.t.len() {
            let zi = (self.t[i] - b) / c;
            let ei = (-0.5 * zi.powi(2)).exp();
            jac[(i, 0)] = -ei;
            jac[(i, 1)] = -(a / c) * ei * zi;
            jac[(i, 2)] = -(a / c) * ei * zi.powi(2);
        }
        Ok(())
    }

    fn fvv(
        &self,
        x: &VecF64,
        v: &VecF64,
        fvv: &mut VecF64,
    ) -> Result<(), Error> {
        let (a, b, c) = (x[0], x[1], x[2]);
        let (va, vb, vc) = (v[0], v[1], v[2]);
        for i in 0..self.t.len() {
            let zi = (self.t[i] - b) / c;
            let ei = (-0.5 * zi.powi(2)).exp();
            let dab = -zi * ei / c;
            let dac = -zi.powi(2) * ei / c;
            let dbb = a * ei / c.powi(2) * (1. - zi.powi(2));
            let dbc = a * zi * ei / c.powi(2) * (2. - zi.powi(2));
            let dcc = a * zi.powi(2) * ei / c.powi(2) * (3. - zi.powi(2));
            fvv[i] = 2.0 * va * vb * dab +
                     2.0 * va * vc * dac +
                           vb * vb * dbb +
                     2.0 * vb * vc * dbc +
                           vc * vc * dcc;
        }
        Ok(())
    }
}

fn solve_system(
    x: &VecF64,
    n: usize,
    fdf: impl Fdf<VecF64>,
    params: &Parameters,
) -> Result<VecF64, Error> {
    let mut work = Workspace::trust(&params, n, x.len())
        .init(x, fdf)?;
    let χsq0 = {
        let f = work.residual();
        dot(&f, &f)?
    };
    // Iterate until convergence.
    println!("            a      b      c    |a|/|v|  cond(J)  |f(x)|");
    work.driver(200)
        .xtol(1e-8).gtol(1e-8).ftol(1e-8)
        .cb(|w, iter| {
            let f = w.residual();
            let x = w.position();
            let avratio = w.avratio();
            let rcond = w.rcond().unwrap();
            println!("iter {:>2}: {:.4} {:.4} {:.4}   {:.4} {:8.4}  {:7.4}",
                iter, x[0], x[1], x[2], avratio, 1.0 / rcond, nrm2(&f));
        })
        .run()?;
    let χsq = {
        let f = work.residual();
        dot(&f, &f)?
    };
    let rcond = work.rcond()?;
    println!("NITER         = {}", work.niter());
    println!("NFEV          = {}", work.nevalf());
    println!("NJEV          = {}", work.nevaldf());
    println!("NAEV          = {}", work.nevalfvv());
    println!("initial cost  = {:.12e}", χsq0);
    println!("final cost    = {:.12e}", χsq);
    println!("final x       = {:?}", work.position());
    println!("final cond(J) = {:.12e}", 1. / rcond);
    Ok(work.position().clone())
}

fn main() -> Result<(), Error> {
    let n = 300;  // number of data points to fit
    let mut params = Parameters::default();
    let (a, b, c) = (5., 0.4, 0.15); // (amplitude, center, width)

    let d = Data::new(a, b, c, n);
    // Starting point
    let x0: VecF64 = [1., 0., 1.].into();
    params.trs = TRS::LMaccel;
    let _x = solve_system(&x0, n, (
        |x: &_, fx: &mut _| d.f(x,fx),
        |x: &_, j: &mut _| d.df(x,j),
        |x: &_, v: &_, fvv: &mut _| d.fvv(x,v,fvv)
        ), &params)?;

    # //assert!(false); // Uncomment to see output
    Ok(())
}
```

### Comparing TRS Methods Example

The following program compares all available nonlinear least squares
trust-region subproblem (TRS) methods on the Branin function, a common
optimization test problem.  The cost function is

$$\begin{aligned}
Φ(x) &= \frac{1}{2} (f₁² + f₂²) \cr
f₁ &= x₂ + a₁ x₁² + a₂ x₁ + a₃ \cr
f₂ &= √{a₄} √{1 + (1 - a₅) \cos x₁}
\end{aligned}$$

with $a₁ = -5.1/(4 π²)$, $a₂ = 5/π$, $a₃ = -6$, $a₄ = 10$, $a₅ =
1/(8π)$.  There are three minima of this function in the range
$(x₁,x₂) ∈ \[-5,15\] × \[-5,15\]$.  The program below uses the
starting point $(x₁, x₂) = (6, 14.5)$ and calculates the solution with
all available nonlinear least squares TRS methods. The program output
is shown below:

```text
Method    NITER NFEV NJEV Initial Cost Final cost Final cond(J) Final x
LM          20   27   21    198.7436     0.39789    6.13992e7    (-3.1416, 12.2750)
LMaccel     27   36   28    198.7436     0.39789    1.44651e7    ( 3.1416,  2.2750)
Dogleg      23   64   23    198.7436     0.39789    5.06919e8    ( 3.1416,  2.2750)
DDogleg     24   69   24    198.7436     0.39789    3.48794e7    ( 3.1416,  2.2750)
Subspace2D  23   54   24    198.7436     0.39789    1.54083e7    ( 3.1416,  2.2750)
```

```text
Method                    NITER  NFEV  NJEV  Initial Cost  Final cost   Final cond(J) Final x
levenberg-marquardt       20     27    21    1.9874e+02    3.9789e-01   6.1399e+07    (-3.14e+00, 1.23e+01)
levenberg-marquardt+accel 27     36    28    1.9874e+02    3.9789e-01   1.4465e+07    (3.14e+00, 2.27e+00)
dogleg                    23     64    23    1.9874e+02    3.9789e-01   5.0692e+08    (3.14e+00, 2.28e+00)
double-dogleg             24     69    24    1.9874e+02    3.9789e-01   3.4879e+07    (3.14e+00, 2.27e+00)
2D-subspace               23     54    24    1.9874e+02    3.9789e-01   2.5142e+07    (3.14e+00, 2.27e+00)
```

The first row of output above corresponds to standard
Levenberg-Marquardt, while the second row includes geodesic
acceleration.  We see that the standard LM method converges to the
minimum at $(-π, 12.275)$ and also uses the least number of iterations
and Jacobian evaluations.  All other methods converge to the minimum
$(π, 2.275)$ and perform similarly in terms of number of Jacobian
evaluations.  We see that $J$ is fairly ill-conditioned at both
minima, indicating that the QR (or SVD) solver is the best choice for
this problem.  Since there are only two parameters in this
optimization problem, we can easily visualize the paths taken by each
method, which are shown below.  The figure shows contours of the cost
function $Φ(x₁, x₂)$ which exhibits three global minima in the range
$\[-5,15\] × \[-5,15\]$.  The paths taken by each solver are shown as
colored lines.

![Paths taken for different TRS methods for the Branin function](https://www.gnu.org/software/gsl/doc/html/_images/nlfit3.png)

The program is given below.

```
use rgsl::{
    Error, MatF64, VecF64,
    blas::d::dot,
    multifit::{Fdf, TRS, Parameters, Workspace},
};
use std::f64::consts::PI;

struct Model {
    a1: f64,
    a2: f64,
    a3: f64,
    a4: f64,
    a5: f64,
}

impl Model {
    fn branin(&self, x: &VecF64, f: &mut VecF64) -> Result<(), Error> {
        let (x1, x2) = (x[0], x[1]);
        f[0] = x2 + self.a1 * x1.powi(2) + self.a2 * x1 + self.a3;
        f[1] = self.a4.sqrt() * (1. + (1. - self.a5) * x1.cos()).sqrt();
        Ok(())
    }

    fn dbranin(&self, x: &VecF64, jac: &mut MatF64) -> Result<(), Error> {
        let x1 = x[0];
        let f2 = self.a4.sqrt() * (1. + (1. - self.a5) * x1.cos()).sqrt();
        jac[(0, 0)] = 2. * self.a1 * x1 + self.a2;
        jac[(0, 1)] = 1.;
        jac[(1, 0)] = -0.5 * self.a4 / f2 * (1. - self.a5) * x1.sin();
        jac[(1, 1)] = 0.;
        Ok(())
    }

    fn branin_vv(&self, x: &VecF64, v: &VecF64, fvv: &mut VecF64) -> Result<(), Error> {
        let x1 = x[0];
        let v1 = v[0];
        let (s, c) = x1.sin_cos();
        let f2 = self.a4.sqrt() * (1. + (1. - self.a5) * c).sqrt();
        let t = 0.5 * self.a4 * (1. - self.a5) / f2;
        fvv[0] = 2. * self.a1 * v1.powi(2);
        fvv[1] = -t * (c + s.powi(2) / f2) * v1.powi(2);
        Ok(())
    }
}

fn solve_system(
    x0: &VecF64,
    params: &Parameters,
    fdf: impl Fdf<VecF64>,
) -> Result<(), Error> {
    let n = 2;
    let work = Workspace::trust(&params, n, x0.len());
    eprintln!("# {:?}/{:?}", work.name(), work.trs_name());
    let mut work = work.init(x0, fdf)?;
    let χsq0 = {
        let f = work.residual();
        dot(&f, &f)?
    };
    work.driver(200)
        .xtol(1e-8).gtol(1e-8).ftol(1e-8)
        .cb(|w, iter| {
            let x = w.position();
            eprintln!("iter {:>2}: {} {}", iter, x[0], x[1])
        })
        .run()?;
    // Store final cost.
    let χsq = {
        let f = work.residual();
        dot(&f, &f)?
    };
    let rcond = work.rcond()?;
    let meth = format!("{:?}", work.trs_name());
    let x = work.position();
    println!("{:11} {:<4} {:<4} {:<5} {:<12.4} {:<10.5} {:<12.5e} ({:>7.4}, {:>7.4})",
        meth, work.niter(),
        work.nevalf(), work.nevaldf(),
        χsq0, χsq,
        1. / rcond,
        x[0], x[1]);
    Ok(())
}

fn main() -> Result<(), Error> {
    let model = Model {
        a1: -5.1 / (4. * PI.powi(2)),
        a2: 5. / PI,
        a3: -6.,
        a4: 10.,
        a5: 1. / (8. * PI),
    };
    let x: VecF64 = [6., 14.5].into();
    println!("{:9} {:<5} {:<4} {:<4} {} {} {} {}",
        "Method", "NITER", "NFEV", "NJEV", "Initial Cost",
        "Final cost", "Final cond(J)", "Final x");
    let mut params = Parameters::default();
    for trs in [TRS::LM, TRS::LMaccel, TRS::Dogleg, TRS::DDogleg, TRS::Subspace2D] {
        params.trs = trs;
        solve_system(&x, &params, (
            |x: &_, fx: &mut _| model.branin(x, fx),
            |x: &_, jac: &mut _| model.dbranin(x, jac),
            |x: &_, v: &_, fvv: &mut _| model.branin_vv(x, v, fvv),
            ))?;
    }
    # //assert!(false); // Uncomment to see output
    Ok(())
}
```

 */

use crate::{
    Error,
    ffi::FFI,
    matrix::{AsMatrix, Matrix, MatrixMut, matrix_as_gsl, matrix_as_gsl_mut},
    vector::{Vector, VectorMut, VectorSpace, vector_as_gsl, vector_as_gsl_mut},
    view::View,
};
use pastey::paste;
use std::{ffi::c_void, ops::ControlFlow, ptr};

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
    /// by using a modified Cholesky decomposition of the matrix
    /// $J^T J + μ D^T D$.  This is more suitable for the dogleg methods
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
    fn as_c(&self) -> sys::gsl_multifit_nlinear_parameters {
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
pub trait Fdf<V: VectorSpace + ?Sized> {
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
    /// [`std::unimplemented`] and let [`Self::has_df`] return `false`.
    /// In this case the Jacobian will be internally computed using
    /// finite difference approximations of the function f.
    fn df(&mut self, x: &V, J: &mut V::Mat) -> Result<(), Error>;

    /// Return `true` iff the method [`Self::df`] is implemented.
    fn has_df(&self) -> bool;

    /// `fvv(x, v, fvv)` must store the `n` components of the vector
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
    V: VectorSpace + ?Sized,
    F: FnMut(&V, &mut V) -> Result<(), Error>,
{
    fn f(&mut self, x: &V, fx: &mut V) -> Result<(), Error> {
        self(x, fx)
    }

    #[inline]
    fn has_df(&self) -> bool {
        false
    }

    fn df(&mut self, _x: &V, _J: &mut V::Mat) -> Result<(), Error> {
        unimplemented!()
    }

    #[inline]
    fn has_fvv(&self) -> bool {
        false
    }

    fn fvv(&mut self, _x: &V, _v: &V, _fvv: &mut V) -> Result<(), Error> {
        unimplemented!()
    }
}

impl<V, F, DF> Fdf<V> for (F, DF)
where
    V: VectorSpace + ?Sized,
    F: FnMut(&V, &mut V) -> Result<(), Error>,
    DF: FnMut(&V, &mut V::Mat) -> Result<(), Error>,
{
    fn f(&mut self, x: &V, fx: &mut V) -> Result<(), Error> {
        self.0(x, fx)
    }

    #[inline]
    fn has_df(&self) -> bool {
        true
    }

    fn df(&mut self, x: &V, J: &mut V::Mat) -> Result<(), Error> {
        self.1(x, J)
    }

    #[inline]
    fn has_fvv(&self) -> bool {
        false
    }

    fn fvv(&mut self, _x: &V, _v: &V, _fvv: &mut V) -> Result<(), Error> {
        unimplemented!()
    }
}

impl<V, F, DF, FVV> Fdf<V> for (F, DF, FVV)
where
    V: VectorSpace + ?Sized,
    F: FnMut(&V, &mut V) -> Result<(), Error>,
    DF: FnMut(&V, &mut V::Mat) -> Result<(), Error>,
    FVV: FnMut(&V, &V, &mut V) -> Result<(), Error>,
{
    fn f(&mut self, x: &V, fx: &mut V) -> Result<(), Error> {
        self.0(x, fx)
    }

    #[inline]
    fn has_df(&self) -> bool {
        true
    }

    fn df(&mut self, x: &V, J: &mut V::Mat) -> Result<(), Error> {
        self.1(x, J)
    }

    #[inline]
    fn has_fvv(&self) -> bool {
        true
    }

    fn fvv(&mut self, x: &V, v: &V, fvv: &mut V) -> Result<(), Error> {
        self.2(x, v, fvv)
    }
}

impl<V, F, DF, FVV> Fdf<V> for (F, Option<DF>, Option<FVV>)
where
    V: VectorSpace + ?Sized,
    F: FnMut(&V, &mut V) -> Result<(), Error>,
    DF: FnMut(&V, &mut V::Mat) -> Result<(), Error>,
    FVV: FnMut(&V, &V, &mut V) -> Result<(), Error>,
{
    fn f(&mut self, x: &V, fx: &mut V) -> Result<(), Error> {
        self.0(x, fx)
    }

    #[inline]
    fn has_df(&self) -> bool {
        self.1.is_some()
    }

    fn df(&mut self, x: &V, J: &mut V::Mat) -> Result<(), Error> {
        self.1.as_mut().expect("df if used")(x, J)
    }

    #[inline]
    fn has_fvv(&self) -> bool {
        self.2.is_some()
    }

    fn fvv(&mut self, x: &V, v: &V, fvv: &mut V) -> Result<(), Error> {
        self.2.as_mut().expect("fvv if used")(x, v, fvv)
    }
}

/// Termination criteria for convergence test.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Test {
    /// The criteria for the relative step size succeeded.
    /// See [`Workspace::test`] for further information.
    Step,
    /// The criteria for the gradient size succeeded.
    /// See [`Workspace::test`] for further information.
    Gradient,
}

impl Test {
    fn from_c(info: i32) -> Test {
        match info {
            1 => Self::Step,
            2 => Self::Gradient,
            _ => panic!(
                "rgsl::multifit::Test: wrong info = {info}, \n\
                         please report to the create owner."
            ),
        }
    }
}

ffi_wrapper!(
    /// Derivative solver.
    #[must_use]
    Workspace<'a, V: VectorSpace + ?Sized>,
    *mut sys::gsl_multifit_nlinear_workspace,
    gsl_multifit_nlinear_free
    ;fdf_struct: Option<Box<sys::gsl_multifit_nlinear_fdf>> => None;
    ;fdf: Option<Box<dyn Fdf<V> + 'a>> => None;
);

ffi_wrapper!(
    /// Uninitialized derivative solver.  Call [`Self::init`] or
    /// [`Self::winit`] to initialize it.
    #[must_use]
    UninitializedWorkspace,
    *mut sys::gsl_multifit_nlinear_workspace,
    gsl_multifit_nlinear_free
);

macro_rules! workspace_common {
    ($fit: ident) => {
        paste! {
            /// Return the number of observations the workspace is for.
            pub fn n(&self) -> usize {
                unsafe {
                    // `f` is non-null because it is allocated by `new`.
                    let w = &* self.unwrap_shared();
                    (*w.f).size
                }
            }

            /// Return the number of parameters the workspace is for.
            pub fn p(&self) -> usize {
                unsafe {
                    // `x` is non-null because it is allocated by `new`.
                    let w = &* self.unwrap_shared();
                    (*w.x).size
                }
            }

            /// Return the name of the solver.
            #[doc(alias = gsl_multi $fit _nlinear_name)]
            pub fn name(&self) -> Type {
                // At the moment, no need to call gsl_multifit_nlinear_name
                // or gsl_multilarge_nlinear_name.  It is part of the type.
                Type::Trust
            }

            /// Return the name of the Trust Region Algorithm used.
            #[doc(alias = gsl_multi $fit _nlinear_trs_name)]
            pub fn trs_name(&self) -> TRS {
                // `trs_name_of_ptr` is implemented differently for
                // `fit` and `large`.
                trs_name_of_ptr(self.unwrap_shared())
            }
        }
    };
}

macro_rules! impl_workspace {
    ($fit: ident, $path: literal) => {
        paste! {
            impl UninitializedWorkspace {
                workspace_common!($fit);

                fn into_workspace<'a, V>(mut self) -> Workspace<'a, V>
                where V: VectorSpace + ?Sized {
                    let work = Workspace::wrap(self.unwrap_unique());
                    // Don't run the destructor on self: we want to
                    // keep the allocated workspace it points to.
                    std::mem::forget(self);
                    work
                }

                /// Initialize the workspace to use the function `fdf`
                /// and the initial guess `x`.
                ///
                /// See [`Workspace::init`] for more information.
                pub fn init<'a, V, F>(
                    self,
                    x: &V,
                    fdf: F,
                ) -> Result<Workspace<'a, V>, Error>
                where
                    V: VectorSpace + ?Sized,
                    F: Fdf<V> + 'a,
                {
                    let mut work = self.into_workspace();
                    work.init(x, fdf)?;
                    Ok(work)
                }

                /// Same as [`Self::init`] but you must in addition
                /// specify a weight vector `w`.  The weighting matrix
                /// is $W = \diag(w₁,w₂,...,wₙ)$.
                pub fn winit<'a, V, F>(
                    self,
                    x: &V,
                    w: &V,
                    fdf: F,
                ) -> Result<Workspace<'a, V>, Error>
                where
                    V: VectorSpace + ?Sized,
                    F: Fdf<V> + 'a,
                {
                    let mut work = self.into_workspace();
                    work.winit(x, w, fdf)?;
                    Ok(work)
                }
            }

            // The associated functions declared in this block do not
            // depend on `V` but Rust needs to know the type of `V`.
            impl<'a> Workspace<'a, crate::vector::VecF64> {
                fn type_to_c(t: Type) -> *const sys::[<gsl_multi $fit _nlinear_type>] {
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
                pub fn new(t: Type, params: &Parameters, n: usize, p: usize) -> UninitializedWorkspace {
                    let w = unsafe { sys::[<gsl_multi $fit _nlinear_alloc>](
                        Self::type_to_c(t), &params.as_c(), n, p)
                    };
                    if w.is_null() {
                        panic!("rgsl::multi{}::Workspace::new: out of memory",
                            stringify!($fit));
                    }
                    UninitializedWorkspace::wrap(w)
                }

                /// Return a new instance of a derivative solver of
                /// type [`Trust`](Type::Trust) for `n` observations
                /// and `p` parameters.
                pub fn trust(params: &Parameters, n: usize, p: usize) -> UninitializedWorkspace {
                    Self::new(Type::Trust, params, n, p)
                }
            }

            impl<'a, V: VectorSpace + ?Sized> Workspace<'a, V> {
                workspace_common!($fit);

                /// Perform a single iteration of the solver.  If the
                /// iteration encounters an unexpected problem then an
                /// error code will be returned.  The solver workspace
                /// maintains a current estimate of the best-fit
                /// parameters at all times.
                ///
                /// See also [`Self::driver`].
                pub fn iterate(&mut self) -> Result<(), Error> {
                    let ret = unsafe {
                        sys::[<gsl_multi $fit _nlinear_iterate>](
                            self.unwrap_unique())
                    };
                    Error::handle(ret, ())
                }

                /// Return the current position (i.e. best-fit
                /// parameters) of the solver, length `p`.
                // checker:ignore
                #[doc(alias = gsl_multi $fit _nlinear_position)]
                pub fn position(&self) -> V::View<'_> {
                    unsafe {
                        let ptr = sys::[<gsl_multi $fit _nlinear_position>](
                            self.unwrap_shared());
                        V::view_from_ptr(ptr)
                    }
                }

                /// Return the function residual vector at the current
                /// position $f(x)$, length `n`.  For weighted
                /// systems, the residual vector includes the
                /// weighting factor $√W$.
                // checker:ignore
                #[doc(alias = gsl_multi $fit _nlinear_residual)]
                pub fn residual(&self) -> V::View<'_> {
                    unsafe {
                        // It is allocated by `new` but irrelevant before `init`.
                        let ptr = sys::[<gsl_multi $fit _nlinear_residual>](
                            self.unwrap_shared());
                        V::view_from_ptr(ptr)
                    }
                }

                /// Return the difference between the current position
                /// and the previous position, i.e. the last step $δ$,
                /// taken as a vector, length `p`.
                pub fn dx(&self) -> V::View<'_> {
                    unsafe {
                        // It is allocated by `new` but irrelevant before `init`.
                        let w = &* self.unwrap_shared();
                        V::view_from_ptr(w.dx)
                    }
                }

                /// Return the number of iterations performed so far.
                /// The iteration counter is updated on each call to
                /// the [`Self::iterate`] function, and reset to 0 in
                /// the [`Self::init`]functions.
                // checker:ignore
                #[doc(alias = gsl_multi $fit _nlinear_niter)]
                pub fn niter(&self) -> usize {
                    // Available after `new` but irrelevant before `init`.
                    unsafe { sys::[<gsl_multi $fit _nlinear_niter>](
                        self.unwrap_shared())
                    }
                }

                /// Return the number of function evaluations.
                pub fn nevalf(&self) -> usize {
                    unsafe {
                        // `fdf` is not null because `Workspace` is
                        // not accessible before `init` is called.
                        let w = &* self.unwrap_shared();
                        (*w.fdf).nevalf
                    }
                }

                /// Return the number of directional derivatives evaluations.
                pub fn nevalfvv(&self) -> usize {
                    unsafe {
                        // `fdf` is not null because `Workspace` is
                        // not accessible before `init` is called.
                        let w = &* self.unwrap_shared();
                        (*w.fdf).nevalfvv
                    }
                }

                /// Return an estimation of the reciprocal condition
                /// number of the Jacobian matrix at the current
                /// position $x$.  The computed value is only an
                /// estimate to give the user a guideline as to the
                /// conditioning of their particular problem.  Its
                /// calculation is based on which factorization method
                /// is used (Cholesky, QR, or SVD).
                ///
                /// - For the Cholesky solver, the matrix $J^T J$ is
                ///   factored at each iteration.  Therefore this
                ///   function will estimate the 1-norm condition number
                ///   $\rcond^2 = 1/(‖J^T J‖₁ · ‖(J^T J)^{-1}‖₁)$
                ///
                /// - For the QR solver, $J$ is factored as $J = Q R$ at
                ///   each iteration.  For simplicity, this function
                ///   calculates the 1-norm conditioning of only the $R$
                ///   factor, $\rcond = 1 / (‖R‖₁ · ‖R^{-1}‖₁)$.  This
                ///   can be computed efficiently since $R$ is upper
                ///   triangular.
                ///
                /// - For the SVD solver, in order to efficiently solve
                ///   the trust region subproblem, the matrix which is
                ///   factored is $J D^{-1}$, instead of $J$ itself.
                ///   The resulting singular values are used to provide
                ///   the 2-norm reciprocal condition number, as $\rcond
                ///   = σ_{\min} / σ_{\max}$.  Note that when using Moré
                ///   scaling, $D ≠ I$ and the resulting $\rcond$
                ///   estimate may be significantly different from the
                ///   true $\rcond$ of $J$ itself.
                pub fn rcond(&self) -> Result<f64, Error> {
                    let mut rcond = f64::NAN;
                    let ret = unsafe {
                        sys::[<gsl_multi $fit _nlinear_rcond>](
                            &mut rcond, self.unwrap_shared())
                    };
                    Error::handle(ret, rcond)
                }

                /// Return the current ratio $|a| / |v|$ of the
                /// acceleration correction term to the velocity step
                /// term.  The acceleration term is computed only by
                /// the [`TRS::LMaccel`] method, so this ratio will be
                /// zero for other TRS methods.
                // checker:ignore
                #[doc(alias = gsl_multi $fit _nlinear_avratio)]
                pub fn avratio(&self) -> f64 {
                    unsafe { sys::[<gsl_multi $fit _nlinear_avratio>](
                        self.unwrap_shared())
                    }
                }

                /// Test for convergence of the minimization method.
                ///
                /// The following criteria is used:
                /// - Testing for a small step size relative to the
                ///   current parameter vector
                ///
                ///   $$|δᵢ| ≤ \xtol (|xᵢ| + \xtol)$$
                ///
                ///   for each $0 ≤ i < p$.  Each element of the step
                ///   vector $δ$ is tested individually in case the
                ///   different parameters have widely different
                ///   scales.  Adding `xtol` to $|xᵢ|$ helps the test
                ///   avoid breaking down in situations where the true
                ///   solution value $xᵢ = 0$.  If this test succeeds,
                ///   the function returns
                ///   [`ControlFlow`]`::Break(Test::Step)`.
                ///
                ///   A general guideline for selecting the step
                ///   tolerance is to choose $\xtol = 10^{-d}$ where
                ///   $d$ is the number of accurate decimal digits
                ///   desired in the solution $x$.  See Dennis and
                ///   Schnabel for more information.
                /// - Testing for a small gradient ($g = ∇Φ(x) = J^T
                ///   f$) indicating a local function minimum:
                ///
                ///   $$\maxᵢ |gᵢ·\max(xᵢ, 1)| ≤ \gtol · \max(Φ(x), 1)$$
                ///
                ///   This expression tests whether the ratio $(∇Φ)ᵢ
                ///   xᵢ / Φ$ is small.  Testing this scaled gradient
                ///   is a better than $∇Φ$ alone since it is a
                ///   dimensionless quantity and so independent of the
                ///   scale of the problem.  The max arguments help
                ///   ensure the test doesn’t break down in regions
                ///   where $xᵢ$ or $Φ(x)$ are close to 0.  If this
                ///   test succeeds, the function returns
                ///   [`ControlFlow`]`::Break(Test::Gradient)`.
                ///
                /// If none of the tests succeed, the function returns
                /// [`ControlFlow`]`::Continue(())`, indicating further
                /// iterations are required.
                ///
                /// See also [`Self::driver`].
                pub fn test(
                    &self,
                    xtol: f64,
                    gtol: f64,
                    ftol: f64,
                ) -> ControlFlow<Test> {
                    let mut info = 0;
                    unsafe {
                        sys::[<gsl_multi $fit _nlinear_test>](
                            xtol,
                            gtol,
                            ftol,
                            &mut info,
                            self.unwrap_shared())
                    };
                    if info == 0 {
                        ControlFlow::Continue(())
                    } else {
                        ControlFlow::Break(Test::from_c(info))
                    }
                }

                /// Provide a high level wrapper that combines the
                /// iteration and convergence testing for easy use.
                ///
                /// Iterate the nonlinear least squares solver for a
                /// maximum of `maxiter` iterations.  After each
                /// iteration, the system is tested for convergence
                /// with the error tolerances `xtol`, `gtol` and
                /// `ftol` (whose default values are `1e-6`.
                /// Additionally, the user may supply a callback
                /// function `cb` which is called after each iteration
                /// as `cb(workspace, iter)`, so that the user may
                /// save or print relevant quantities for each
                /// iteration.  Upon successful convergence, the
                /// function returns `Ok(test)` with `test` being the
                /// reason for convergence (see [`Self::test`]).
                ///
                /// If the function has not converged after `maxiter`
                /// iterations, `Err(Error::MaxIteration)` is
                /// returned.  (Note that this can happen if the
                /// derivatives of $f$ are incorrect.)  In rare cases,
                /// during an iteration the algorithm may be unable to
                /// find a new acceptable step $δ$ to take.  In this
                /// case, `Err(Error::NoProgress)` is returned
                /// indicating no further progress can be made.  If
                /// your problem is having difficulty converging, see
                /// [Troubleshooting](crate::multifit#Troubleshooting)
                /// for further guidance.
                ///
                /// # Example
                ///
                /// ```
                /// use rgsl::{Error, multifit::{Workspace, Parameters}, VecF64};
                /// let p = Parameters::default();
                /// let x: VecF64 = [10.].into();
                /// let f = |x: &VecF64, fx: &mut VecF64| {
                ///     fx[0] = x[0];
                ///     fx[1] = x[0] - 2.;
                ///     Ok(())
                /// };
                /// let mut w = Workspace::trust(&p, 2, 1)
                ///     .init(&x, f)?;
                /// w.driver(10)
                ///     .cb(|_, iter| println!("iter: {iter}"))
                ///     .run()?;
                /// assert!((w.position()[0] - 1.) < 1e-9);
                /// # Ok::<(), Error>(())
                /// ```
                pub fn driver<'b>(&'b mut self, maxiter: usize) -> Driver<'a, 'b, V> {
                    Driver {
                        w: self,
                        maxiter,
                        xtol: 1e-6,
                        gtol: 1e-6,
                        ftol: 1e-6,
                        cb: None,
                    }
                }
            }

            #[must_use]
            pub struct Driver<'a, 'b, V: VectorSpace + ?Sized> {
                w: &'b mut Workspace<'a, V>,
                maxiter: usize,
                xtol: f64,
                gtol: f64,
                ftol: f64,
                cb: Option<Box<dyn FnMut(&Workspace<'a, V>, usize) + 'a>>,
            }

            impl<'a, 'b, V: VectorSpace + ?Sized> Driver<'a, 'b, V> {
                /// Set the relative tolerance for the step size.
                ///
                /// See [`Workspace::test`] for more information.
                pub fn xtol(mut self, xtol: f64) -> Self {
                    self.xtol = xtol;
                    self
                }

                /// Set the tolerance for the gradient.
                ///
                /// See [`Workspace::test`] for more information.
                pub fn gtol(mut self, gtol: f64) -> Self {
                    self.gtol = gtol;
                    self
                }

                /// Set the tolerance for function value.
                pub fn ftol(mut self, ftol: f64) -> Self {
                    self.ftol = ftol;
                    self
                }

                /// Add a callback `cb` to be executed at each iteration.
                pub fn cb(
                    mut self,
                    cb: impl FnMut(&Workspace<V>, usize) + 'a
                ) -> Self {
                    self.cb = Some(Box::new(cb));
                    self
                }

                /// Run the driver and return the reason for stopping.
                ///
                /// See [`test`] for more information.
                pub fn run(self) -> Result<Test, Error> {
                    type F<'a, 'c, V> = Box<dyn FnMut(&Workspace<'a, V>, usize) + 'c>;
                    type P<'c> = &'c mut bool;

                    unsafe extern "C" fn trampoline_cb<V: VectorSpace + ?Sized>(
                        iter: usize,
                        params: *mut c_void,
                        w: *const sys::[<gsl_multi $fit _nlinear_workspace>],
                    ) {
                        let ret = std::panic::catch_unwind(|| unsafe {
                            let (cb, _) = &mut *params.cast::<(F<V>, P)>();
                            // `View` doesn't allow mutable access so
                            // the cast is fine.
                            let w = Workspace::<V>::wrap(w as *mut _);
                            // We do not own the workspace and C
                            // memory must not be freed when it is
                            // dropped.
                            let w = View::new(w, false);
                            cb(&w, iter)
                        });
                        if ret.is_err() {
                            let panicked = unsafe {
                                &mut (&mut *params.cast::<(F<V>, P)>()).1
                            };
                            **panicked = true;
                        }
                    }

                    let mut info = 0;
                    let ret = if let Some(cb) = self.cb {
                        let mut panicked = false;
                        let mut params = (cb, &mut panicked);
                        let ret = unsafe { sys::[<gsl_multi $fit _nlinear_driver>](
                            self.maxiter, self.xtol, self.gtol, self.ftol,
                            Some(trampoline_cb::<V>),
                            &mut params as *mut (F<'a, 'b, V>, P) as *mut _,
                            &mut info,
                            self.w.unwrap_unique(),
                        )};
                        if panicked {
                            return Err(Error::BadFunction)
                        }
                        ret
                    } else {
                        unsafe { sys::[<gsl_multi $fit _nlinear_driver>](
                            self.maxiter, self.xtol, self.gtol, self.ftol,
                            None,
                            ptr::null_mut(),
                            &mut info,
                            self.w.unwrap_unique(),
                        )}
                    };
                    match Error::handle(ret, ()) {
                        Err(e) => Err(e),
                        // Only evaluate `from_c` if no error.
                        Ok(_) => Ok(Test::from_c(info)),
                    }
                }
            }
        }
    };
}

impl_workspace!(fit, "rgsl::multifit");

// Common to `UninitializedWorkspace` and `Workspace` but differs for
// the "large" version.
#[doc(alias = "gsl_multifit_nlinear_trs_name")]
// checker:ignore
fn trs_name_of_ptr(ptr: *const sys::gsl_multifit_nlinear_workspace) -> TRS {
    let n = unsafe { sys::gsl_multifit_nlinear_trs_name(ptr) };
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

impl<'a, V: VectorSpace + ?Sized> Workspace<'a, V> {
    // checker:ignore
    #[doc(alias = "gsl_multifit_nlinear_winit")]
    fn init_wts<F: Fdf<V> + 'a>(&mut self, x: &V, wts: Option<&V>, fdf: F) -> Result<(), Error> {
        unsafe extern "C" fn f_trampoline<V: VectorSpace + ?Sized, F: Fdf<V>>(
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
            Error::to_c(ret.map_err(|_| Error::BadFunction).flatten())
        }

        unsafe extern "C" fn df_trampoline<V: VectorSpace + ?Sized, F: Fdf<V>>(
            x: *const sys::gsl_vector,
            params: *mut c_void,
            J: *mut sys::gsl_matrix,
        ) -> i32 {
            let ret = std::panic::catch_unwind(|| unsafe {
                let fdf = &mut *params.cast::<F>();
                let vx = V::view_from_ptr(x);
                let mut vJ = V::Mat::mat_view_from_mut_ptr(J);
                fdf.df(&*vx, &mut vJ)
            });
            Error::to_c(ret.map_err(|_| Error::BadFunction).flatten())
        }

        unsafe extern "C" fn fvv_trampoline<V: VectorSpace + ?Sized, F: Fdf<V>>(
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
            Error::to_c(ret.map_err(|_| Error::BadFunction).flatten())
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

        let x = vector_as_gsl(x);
        let wts = if let Some(wts) = wts {
            &*vector_as_gsl(wts)
        } else {
            ptr::null()
        };
        let ret = unsafe {
            sys::gsl_multifit_nlinear_winit(
                &*x,              // copied into the workspace
                wts,              // copied into the workspace
                &mut *fdf_struct, // Heap pointer (stable address)
                self.unwrap_unique(),
            )
        };
        self.fdf_struct = Some(fdf_struct);
        Error::handle(ret, ())
    }

    /// Reinitialize the workspace to use the function `fdf` and the
    /// initial guess `x`.
    ///
    /// There are several ways of specifying a function.
    /// - One can simply pass a function `f` such that the call `f(x,
    ///   fx)` stores the `n` components of the vector $f(x)$ in `fx`
    ///   for argument `x`, returning an appropriate error code if the
    ///   function cannot be computed.  This function must always be
    ///   passed.
    /// - One can pass `(f, df)` where `f` is as above and `df(x, J)`
    ///   stores the `n`-by-`p` Jacobian in `J`, see [`Fdf::df`].
    /// - One can give `(f, df, fvv)` where `f` and `df` are as above
    ///   and `fvv(x, v, fvv)` stores the second directional
    ///   derivative in `fvv`, see [`Fdf::fvv`].
    /// - One can pass `(f, df_opt, fvv_opt)` where `f` is as above,
    ///   `df_opt` is `Some(df)` with `df` as above or `None`, and
    ///   `fvv_opt` is `Some(fvv)` with `fvv` as above or `None`.
    ///
    /// If one of the functions `f`, `df`, or `fvv` panics, the error
    /// [`Error::BadFunction`] is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use rgsl::{Error, multifit::{Workspace, Parameters}, VecF64};
    /// let p = Parameters::default();
    /// let x: VecF64 = [10.].into();
    /// let f = |x: &VecF64, fx: &mut VecF64| {
    ///     fx[0] = x[0];
    ///     fx[1] = x[0] - 2.;
    ///     Ok(())
    /// };
    /// let mut w = Workspace::trust(&p, 2, 1)
    ///     .init(&x, f)?;
    /// while ! matches!(w.iterate(), Err(Error::NoProgress)) {}
    /// assert!((w.position()[0] - 1.) < 1e-8);
    /// let jac = w.jac();
    /// assert!((jac[(0,0)] - 1.) < 1e-15);
    /// assert!((jac[(1,0)] - 1.) < 1e-15);
    /// # Ok::<(), Error>(())
    /// ```
    #[doc(alias = "gsl_multifit_nlinear_init")]
    pub fn init<F: Fdf<V> + 'a>(&mut self, x: &V, fdf: F) -> Result<(), Error> {
        self.init_wts(x, None, fdf)
    }

    /// Same as [`Self::init`] but you must in addition specify a
    /// weight vector `w`.  The weighting matrix is $W =
    /// \diag(w₁,w₂,...,wₙ)$.
    #[doc(alias = "gsl_multifit_nlinear_winit")]
    pub fn winit<F: Fdf<V> + 'a>(&mut self, x: &V, w: &V, fdf: F) -> Result<(), Error> {
        self.init_wts(x, Some(w), fdf)
    }

    /// Return the Jacobian matrix at the current position $J(x)$,
    /// size `n`-by-`p`.
    #[doc(alias = "gsl_multifit_nlinear_jac")]
    pub fn jac(&self) -> <V::Mat as AsMatrix>::MatView<'_> {
        unsafe {
            let ptr = sys::gsl_multifit_nlinear_jac(self.unwrap_shared());
            debug_assert!(!ptr.is_null());
            V::Mat::mat_view_from_ptr(ptr)
        }
    }

    /// Return the number of Jacobian evaluations.
    pub fn nevaldf(&self) -> usize {
        // See `nevalf`.
        unsafe {
            let w = &*self.unwrap_shared();
            (*w.fdf).nevaldf
        }
    }
}

#[cfg(feature = "v2_3")]
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
        fn as_c(&self) -> sys::gsl_multilarge_nlinear_parameters {
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
    pub trait Fdf<V: VectorSpace + ?Sized> {
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
            JTJ: Option<&mut V::Mat>,
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

    impl<V, F, DF> Fdf<V> for (F, DF)
    where
        V: VectorSpace + ?Sized,
        F: FnMut(&V, &mut V) -> Result<(), Error>,
        DF: FnMut(bool, &V, &V, &mut V, Option<&mut V::Mat>) -> Result<(), Error>,
    {
        fn f(&mut self, x: &V, fx: &mut V) -> Result<(), Error> {
            self.0(x, fx)
        }

        fn df(
            &mut self,
            transJ: bool,
            x: &V,
            u: &V,
            v: &mut V,
            JTJ: Option<&mut V::Mat>,
        ) -> Result<(), Error> {
            self.1(transJ, x, u, v, JTJ)
        }

        #[inline]
        fn has_fvv(&self) -> bool {
            false
        }

        fn fvv(&mut self, _x: &V, _v: &V, _fvv: &mut V) -> Result<(), Error> {
            unimplemented!()
        }
    }

    impl<V, F, DF, FVV> Fdf<V> for (F, DF, FVV)
    where
        V: VectorSpace + ?Sized,
        F: FnMut(&V, &mut V) -> Result<(), Error>,
        DF: FnMut(bool, &V, &V, &mut V, Option<&mut V::Mat>) -> Result<(), Error>,
        FVV: FnMut(&V, &V, &mut V) -> Result<(), Error>,
    {
        fn f(&mut self, x: &V, fx: &mut V) -> Result<(), Error> {
            self.0(x, fx)
        }

        fn df(
            &mut self,
            transJ: bool,
            x: &V,
            u: &V,
            v: &mut V,
            JTJ: Option<&mut V::Mat>,
        ) -> Result<(), Error> {
            self.1(transJ, x, u, v, JTJ)
        }

        #[inline]
        fn has_fvv(&self) -> bool {
            true
        }

        fn fvv(&mut self, x: &V, v: &V, fvv: &mut V) -> Result<(), Error> {
            self.2(x, v, fvv)
        }
    }

    impl<V, F, DF, FVV> Fdf<V> for (F, DF, Option<FVV>)
    where
        V: VectorSpace + ?Sized,
        F: FnMut(&V, &mut V) -> Result<(), Error>,
        DF: FnMut(bool, &V, &V, &mut V, Option<&mut V::Mat>) -> Result<(), Error>,
        FVV: FnMut(&V, &V, &mut V) -> Result<(), Error>,
    {
        fn f(&mut self, x: &V, fx: &mut V) -> Result<(), Error> {
            self.0(x, fx)
        }

        fn df(
            &mut self,
            transJ: bool,
            x: &V,
            u: &V,
            v: &mut V,
            JTJ: Option<&mut V::Mat>,
        ) -> Result<(), Error> {
            self.1(transJ, x, u, v, JTJ)
        }

        #[inline]
        fn has_fvv(&self) -> bool {
            self.2.is_some()
        }

        fn fvv(&mut self, x: &V, v: &V, fvv: &mut V) -> Result<(), Error> {
            self.2.as_mut().expect("fvv if used")(x, v, fvv)
        }
    }

    ffi_wrapper!(
        /// Derivative solver for large problems.
        #[must_use]
        Workspace<'a, V: VectorSpace + ?Sized>,
        *mut sys::gsl_multilarge_nlinear_workspace,
        gsl_multilarge_nlinear_free
        ;fdf_struct: Option<Box<sys::gsl_multilarge_nlinear_fdf>> => None;
        ;fdf: Option<Box<dyn Fdf<V> + 'a>> => None;
    );

    ffi_wrapper!(
        /// Uninitialized derivative solver.  Call [`Self::init`] or
        /// [`Self::winit`] to initialize it.
        #[must_use]
        UninitializedWorkspace,
        *mut sys::gsl_multilarge_nlinear_workspace,
        gsl_multilarge_nlinear_free
    );

    fn trs_name_of_ptr(ptr: *const sys::gsl_multilarge_nlinear_workspace) -> TRS {
        let n = unsafe { sys::gsl_multilarge_nlinear_trs_name(ptr) };
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

    impl_workspace!(large, "rgsl::multilarge");

    impl<'a, V: VectorSpace + ?Sized> Workspace<'a, V> {
        #[doc(alias = "gsl_multilarge_nlinear_winit")]
        fn init_wts<F: Fdf<V> + 'a>(
            &mut self,
            x: &V,
            wts: Option<&V>,
            fdf: F,
        ) -> Result<(), Error> {
            unsafe extern "C" fn f_trampoline<V: VectorSpace + ?Sized, F: Fdf<V>>(
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
                Error::to_c(ret.map_err(|_| Error::BadFunction).flatten())
            }

            unsafe extern "C" fn df_trampoline<V: VectorSpace + ?Sized, F: Fdf<V>>(
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
                    if JTJ.is_null() {
                        fdf.df(transJ, &*vx, &*vu, &mut *vv, None)
                    } else {
                        let mut JTJ = V::Mat::mat_view_from_mut_ptr(JTJ);
                        fdf.df(transJ, &*vx, &*vu, &mut *vv, Some(&mut *JTJ))
                    }
                });
                Error::to_c(ret.map_err(|_| Error::BadFunction).flatten())
            }

            unsafe extern "C" fn fvv_trampoline<V: VectorSpace + ?Sized, F: Fdf<V>>(
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
                Error::to_c(ret.map_err(|_| Error::BadFunction).flatten())
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

            let x = vector_as_gsl(x);
            let wts = if let Some(wts) = wts {
                &*vector_as_gsl(wts)
            } else {
                ptr::null()
            };
            let ret = unsafe {
                sys::gsl_multilarge_nlinear_winit(
                    &*x,              // copied into the workspace
                    wts,              // copied into the workspace
                    &mut *fdf_struct, // Heap pointer (stable address)
                    self.unwrap_unique(),
                )
            };
            self.fdf_struct = Some(fdf_struct);
            Error::handle(ret, ())
        }

        /// Initialize, or reinitialize, the workspace to use the function
        /// `fdf` and the initial guess `x`.
        ///
        /// There are several ways of specifying a function.
        /// - One can pass `(f, df)` where `f(x, fx)` stores the `n`
        ///   components of the vector $f(x)$ in `fx` for argument
        ///   `x`, returning an appropriate error code if the function
        ///   cannot be computed (see [`Fdf::f`]) and `df(transJ, x,
        ///   u, v)` stores the product $Ju ∈ ℝⁿ$, where $u ∈ ℝᵖ$, in
        ///   `v`, see [`Fdf::df`].
        /// - One can give `(f, df, fvv)` where `f` and `df` are as above
        ///   and `fvv(x, v, fvv)` stores the second directional
        ///   derivative in `fvv`, see [`Fdf::fvv`].
        /// - One can pass `(f, df_opt, fvv_opt)` where `f` is as above,
        ///   `df_opt` is `Some(df)` with `df` as above or `None`, and
        ///   `fvv_opt` is `Some(fvv)` with `fvv` as above or `None`.
        ///
        /// If one of the functions `f`, `df`, or `fvv` panics, the error
        /// [`Error::BadFunction`] is returned.
        ///
        /// # Examples
        ///
        /// ```
        /// use rgsl::{multilarge::{Workspace, Parameters}, VecF64, MatF64};
        /// let p = Parameters::default();
        /// let mut w = Workspace::trust(&p, 2, 1);
        /// let x = VecF64::from_slice(&[1.]);
        /// let f = |x: &VecF64, fx: &mut VecF64| {
        ///     fx[0] = x[0];
        ///     fx[1] = x[0] - 2.;
        ///     Ok(())
        /// };
        /// let df = |transJ, x: &VecF64, u: &VecF64, v:&mut VecF64, JTJ: Option<&mut MatF64>| {
        ///     if transJ {
        ///         v[0] = u[0] + u[1];
        ///         if let Some(JTJ) = JTJ {
        ///             JTJ[(0,0)] = 2.;
        ///         }
        ///     } else {
        ///         v[0] = u[0];
        ///         v[1] = u[0];
        ///     }
        ///     Ok(())
        /// };
        /// w.init(&x, (f, df));
        /// ```
        #[doc(alias = "gsl_multilarge_nlinear_init")]
        pub fn init<F: Fdf<V> + 'a>(&mut self, x: &V, fdf: F) -> Result<(), Error> {
            self.init_wts(x, None, fdf)
        }

        /// Same as [`Self::init`] except that one can specify a weight
        /// vector `w`.
        #[doc(alias = "gsl_multilarge_nlinear_winit")]
        pub fn winit<F: Fdf<V> + 'a>(&mut self, x: &V, w: &V, fdf: F) -> Result<(), Error> {
            self.init_wts(x, Some(w), fdf)
        }

        /// Return the number of matrix-vector $Ju$ evaluations.
        pub fn nevaldf(&self) -> usize {
            // See `nevalf`.
            unsafe {
                let w = &*self.unwrap_shared();
                (*w.fdf).nevaldfu
            }
        }

        /// Compute the covariance matrix of best-fit parameters and
        /// stores it in `covar`.
        ///
        /// See [`covar`] for more information.
        #[doc(alias = "gsl_multilarge_nlinear_covar")]
        pub fn covar(&mut self, covar: &mut (impl MatrixMut<f64> + ?Sized)) -> Result<(), Error> {
            let mut covar = matrix_as_gsl_mut(covar);
            let ret =
                unsafe { sys::gsl_multilarge_nlinear_covar(&mut *covar, self.unwrap_unique()) };
            Error::handle(ret, ())
        }
    }
}

/// Compute the covariance matrix of best-fit parameters using the
/// Jacobian matrix $J$ and stores it in `covar`.
///
/// The parameter `epsrel` is used to remove linear-dependent columns
/// when $J$ is rank deficient.
///
/// The covariance matrix is given by,
///
/// $$C = (J^T J)^{-1}$$
///
/// or in the weighted case,
///
/// C = (J^T W J)^{-1}
///
/// and is computed using the factored form of the Jacobian (Cholesky,
/// QR, or SVD).  Any columns of $R$ which satisfy
///
/// $$|R_{kk}| ≤ \epsrel |R_{11}|$$
///
/// are considered linearly-dependent and are excluded from the
/// covariance matrix (the corresponding rows and columns of the
/// covariance matrix are set to zero).
///
/// If the minimisation uses the weighted least-squares function $fᵢ =
/// (Y(x, tᵢ) - yᵢ) / σᵢ$ then the covariance matrix above gives the
/// statistical error on the best-fit parameters resulting from the
/// Gaussian errors $σᵢ$ on the underlying data $yᵢ$.  This can be
/// verified from the relation $δf = J δc$ and the fact that the
/// fluctuations in $f$ from the data $yᵢ$ are normalised by $σᵢ$ and
/// so satisfy
///
/// $$⟨δf δf^T⟩ = I.$$
///
/// For an unweighted least-squares function $fᵢ = (Y(x, tᵢ) - yᵢ)$
/// the covariance matrix above should be multiplied by the variance
/// of the residuals about the best-fit $σ² = ∑ (yᵢ - Y(x,tᵢ))² /
/// (n-p)$ to give the variance-covariance matrix $σ² C$.  This
/// estimates the statistical error on the best-fit parameters from
/// the scatter of the underlying data.
///
/// For more information about covariance matrices see [Linear
/// Least-Squares Overview](crate::fit#overview).
#[doc(alias = "gsl_multifit_nlinear_covar")]
pub fn covar(
    J: &(impl Matrix<f64> + ?Sized),
    epsrel: f64,
    covar: &mut (impl MatrixMut<f64> + ?Sized),
) -> Result<(), Error> {
    let J = matrix_as_gsl(J);
    let mut covar = matrix_as_gsl_mut(covar);
    let ret = unsafe { sys::gsl_multifit_nlinear_covar(&*J, epsrel, &mut *covar) };
    Error::handle(ret, ())
}

#[doc(alias = "gsl_multifit_test_delta")]
pub fn test_delta(
    dx: &(impl Vector<f64> + ?Sized),
    x: &(impl Vector<f64> + ?Sized),
    epsabs: f64,
    epsrel: f64,
) -> Result<(), Error> {
    let dx = vector_as_gsl(dx);
    let x = vector_as_gsl(x);
    let ret = unsafe { sys::gsl_multifit_test_delta(&*dx, &*x, epsabs, epsrel) };
    Error::handle(ret, ())
}

#[doc(alias = "gsl_multifit_gradient")]
pub fn gradient(
    J: &(impl Matrix<f64> + ?Sized),
    f: &(impl Vector<f64> + ?Sized),
    g: &mut (impl VectorMut<f64> + ?Sized),
) -> Result<(), Error> {
    let J = matrix_as_gsl(J);
    let f = vector_as_gsl(f);
    let mut g = vector_as_gsl_mut(g);
    let ret = unsafe { sys::gsl_multifit_gradient(&*J, &*f, &mut *g) };
    Error::handle(ret, ())
}

#[doc(alias = "gsl_multifit_linear_lreg")]
pub fn linear_lreg(
    smin: f64,
    smax: f64,
    reg_param: &mut (impl VectorMut<f64> + ?Sized),
) -> Result<(), Error> {
    let mut reg_param = vector_as_gsl_mut(reg_param);
    let ret = unsafe { sys::gsl_multifit_linear_lreg(smin, smax, &mut *reg_param) };
    Error::handle(ret, ())
}

/// Returns `idx`.
#[doc(alias = "gsl_multifit_linear_lcorner")]
pub fn linear_lcorner(
    rho: &(impl Vector<f64> + ?Sized),
    eta: &(impl Vector<f64> + ?Sized),
) -> Result<usize, Error> {
    let mut idx = 0;
    let rho = vector_as_gsl(rho);
    let eta = vector_as_gsl(eta);
    let ret = unsafe { sys::gsl_multifit_linear_lcorner(&*rho, &*eta, &mut idx) };
    Error::handle(ret, idx)
}

/// Returns `(Value, idx)`.
#[doc(alias = "gsl_multifit_linear_lcorner2")]
pub fn linear_lcorner2(
    rho: &(impl Vector<f64> + ?Sized),
    eta: &(impl Vector<f64> + ?Sized),
) -> Result<usize, Error> {
    let mut idx = 0;
    let rho = vector_as_gsl(rho);
    let eta = vector_as_gsl(eta);
    let ret = unsafe { sys::gsl_multifit_linear_lcorner2(&*rho, &*eta, &mut idx) };
    Error::handle(ret, idx)
}

#[doc(alias = "gsl_multifit_linear_Lk")]
pub fn linear_Lk(p: usize, k: usize, L: &mut (impl MatrixMut<f64> + ?Sized)) -> Result<(), Error> {
    let mut L = matrix_as_gsl_mut(L);
    let ret = unsafe { sys::gsl_multifit_linear_Lk(p, k, &mut *L) };
    Error::handle(ret, ())
}

#[doc(alias = "gsl_multilarge_linear_L_decomp")]
pub fn linear_L_decomp(
    L: &mut (impl MatrixMut<f64> + ?Sized),
    tau: &mut (impl VectorMut<f64> + ?Sized),
) -> Result<(), Error> {
    let mut L = matrix_as_gsl_mut(L);
    let mut tau = vector_as_gsl_mut(tau);
    let ret = unsafe { sys::gsl_multilarge_linear_L_decomp(&mut *L, &mut *tau) };
    Error::handle(ret, ())
}

#[cfg(any(test, doctest))]
mod test {
    /// One cannot use `iterate` if the workspace was not initialized.
    ///
    /// ```compile_fail
    /// use rgsl::{Error, multifit::{Parameters, Workspace}};
    /// let params = Parameters::default();
    /// let mut w = Workspace::trust(&params, 1, 1);
    /// w.iterate();
    /// ```
    ///
    /// Also, the result of `w.init` cannot be ignored.
    ///
    /// ```compile_fail
    /// #![deny(unused_must_use)]
    /// use rgsl::{Error, multifit::{Parameters, Workspace}, VecF64};
    /// let params = Parameters::default();
    /// let x0: VecF64 = [1.].into();
    /// let w = Workspace::trust(&params, 1, x0.len());
    /// // Return value must not be ignored.
    /// w.init(&x0, |x: &VecF64, fx: &mut VecF64| {
    ///     fx[0] = x[0];
    ///     Ok(())
    /// })?;
    /// Ok::<(), Error>(())
    /// ```
    use super::*;
    use crate::VecF64;

    #[test]
    fn test_name() {
        type W<'a> = Workspace<'a, VecF64>;

        let params = Parameters::default();
        let w = W::trust(&params, 10, 2);
        assert_eq!(w.trs_name(), TRS::LM);

        let p = Parameters {
            trs: TRS::LMaccel,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).trs_name(), TRS::LMaccel);

        let p = Parameters {
            trs: TRS::Dogleg,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).trs_name(), TRS::Dogleg);

        let p = Parameters {
            trs: TRS::DDogleg,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).trs_name(), TRS::DDogleg);

        let p = Parameters {
            trs: TRS::Subspace2D,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).trs_name(), TRS::Subspace2D);
    }

    #[test]
    fn test_large_name() {
        use large::{Parameters, TRS};
        type W<'a> = large::Workspace<'a, VecF64>;

        let params = Parameters::default();
        let w = W::trust(&params, 10, 2);
        assert_eq!(w.trs_name(), TRS::LM);

        let p = Parameters {
            trs: TRS::LMaccel,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).trs_name(), TRS::LMaccel);

        let p = Parameters {
            trs: TRS::Dogleg,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).trs_name(), TRS::Dogleg);

        let p = Parameters {
            trs: TRS::DDogleg,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).trs_name(), TRS::DDogleg);

        let p = Parameters {
            trs: TRS::Subspace2D,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).trs_name(), TRS::Subspace2D);

        let p = Parameters {
            trs: TRS::CgST,
            ..params
        };
        assert_eq!(W::trust(&p, 10, 2).trs_name(), TRS::CgST);
    }

    #[test]
    fn panicking_fn() -> Result<(), Error> {
        let params = Parameters::default();
        let w = Workspace::trust(&params, 1, 1);
        let x0: VecF64 = [1.].into();
        let ret = w.init(&x0, |_x: &VecF64, _fx: &mut VecF64| {
            panic!("The function `f` panics.");
        });
        assert!(matches!(ret, Err(Error::BadFunction)));

        Ok(())
    }
}
