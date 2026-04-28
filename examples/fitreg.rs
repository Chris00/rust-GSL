//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

#[cfg(feature = "v2_5")]
mod example {
    use rgsl::{MatF64, MultifitLinearWorkspace, Rng, RngType, VecF64, multifit};

    const N: usize = 1000; // number of observations
    const P: usize = 2; // number of model parameters
    const NPOINTS: usize = 200; // number of points on L-curve and GCV curve

    pub fn run() {
        let mut r = Rng::new(RngType::default()).expect("Rng::new failed");
        let mut x = MatF64::new(N, P);
        let mut y = VecF64::new(N);

        for i in 0..N {
            // generate first random variable u
            let ui = 5. * r.gaussian(1.);

            // set v = u + noise
            let vi = ui + r.gaussian(0.001);

            // set y = u + v + noise
            let yi = ui + vi + r.gaussian(1.);

            // since u =~ v, the matrix X is ill-conditioned
            x.set(i, 0, ui);
            x.set(i, 1, vi);

            // rhs vector
            y.set(i, yi);
        }

        let mut w =
            MultifitLinearWorkspace::new(N, P).expect("MultifitLinearWorkspace::new failed");
        // OLS solution
        let mut c = VecF64::new(P);
        // regularized solution (L-curve)
        let mut c_lcurve = VecF64::new(P);
        // regularized solution (GCV)
        let mut c_gcv = VecF64::new(P);

        let mut reg_param = VecF64::new(NPOINTS);
        // residual norms
        let mut rho = VecF64::new(NPOINTS);
        // solution norms
        let mut eta = VecF64::new(NPOINTS);
        // GCV function values
        let mut g = VecF64::new(NPOINTS);

        // compute SVD of X
        w.linear_svd(&mut x).unwrap();

        // Get reciprocal condition number of X
        let rcond = w.linear_rcond();
        eprintln!("matrix condition number = {}", 1. / rcond);
        eprintln!();

        // unregularized (standard) least squares fit, lambda = 0
        let (rnorm, snorm) = w.linear_solve(0., &x, &y, &mut c).unwrap();
        let chisq = rnorm.powi(2);

        eprintln!("\n=== Unregularized fit ===");
        eprintln!("best fit: y = {} u + {} v", c.get(0), c.get(1));
        eprintln!("residual norm = {}", rnorm);
        eprintln!("solution norm = {}", snorm);
        eprintln!("chisq/dof = {}", chisq / (N - P) as f64);

        // calculate L-curve and find its corner
        w.linear_lcurve(&y, &mut reg_param, &mut rho, &mut eta)
            .unwrap();
        let reg_idx = multifit::linear_lcorner(&rho, &eta).unwrap();

        // store optimal regularization parameter
        let lambda_l = reg_param.get(reg_idx);

        // regularize with lambda_l
        let (rnorm, snorm) = w.linear_solve(lambda_l, &x, &y, &mut c_lcurve).unwrap();
        let chisq = rnorm.powi(2) + (lambda_l * snorm).powi(2);

        eprintln!("\n=== Regularized fit (L-curve) ===");
        eprintln!("optimal lambda: {}", lambda_l);
        eprintln!(
            "best fit: y = {} u + {} v",
            c_lcurve.get(0),
            c_lcurve.get(1)
        );
        eprintln!("residual norm = {}", rnorm);
        eprintln!("solution norm = {}", snorm);
        eprintln!("chisq/dof = {}", chisq / (N - P) as f64);

        // calculate GCV curve and find its minimum
        let (lambda_gcv, g_gcv) = w.linear_gcv(&y, &mut reg_param, &mut g).unwrap();

        // regularize with lambda_gcv
        let (rnorm, snorm) = w.linear_solve(lambda_gcv, &x, &y, &mut c_gcv).unwrap();
        let chisq = rnorm.powi(2) + (lambda_gcv * snorm).powi(2);

        eprintln!("\n=== Regularized fit (GCV) ===\n");
        eprintln!("optimal lambda: {}", lambda_gcv);
        eprintln!("best fit: y = {} u + {} v", c_gcv.get(0), c_gcv.get(1));
        eprintln!("residual norm = {}", rnorm);
        eprintln!("solution norm = {}", snorm);
        eprintln!("chisq/dof = {}", chisq / (N - P) as f64);

        // output L-curve and GCV curve
        for i in 0..NPOINTS {
            println!(
                "{:.6} {:.6} {:.6} {:.6}",
                reg_param.get(i),
                rho.get(i),
                eta.get(i),
                g.get(i)
            );
        }

        // output L-curve corner point
        println!("\n\n{:.6} {:.6}", rho.get(reg_idx), eta.get(reg_idx));

        // output GCV curve corner minimum
        println!("\n\n{:.6} {:.6}", lambda_gcv, g_gcv);
    }
}

#[cfg(feature = "v2_5")]
fn main() {
    example::run();
}

#[cfg(not(feature = "v2_5"))]
fn main() {
    eprintln!("You need to enable the `v2_5` feature to be able to run this example");
}
