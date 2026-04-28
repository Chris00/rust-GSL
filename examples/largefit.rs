//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#[cfg(feature = "v2_2")]
mod example {
    use rgsl::{
        MatF64, MultilargeLinearType, MultilargeLinearWorkspace, Rng, RngType, VecF64, blas,
    };

    // number of observations
    const N: usize = 50000;
    // polynomial order + 1
    const P: usize = 16;
    // regularization parameter
    const LAMBDA: f64 = 0.;
    // number of blocks to accumulate
    const NBLOCK: usize = 5;
    // number of rows per block
    const NROWS: usize = N / NBLOCK;
    const NLCURVE: usize = 200;
    const DT: f64 = 1. / (N as f64 - 1.);

    fn func(t: f64) -> f64 {
        let x = (10. * t).sin();
        (x * x * x).exp()
    }

    fn build_row(t: f64, row: &mut VecF64) {
        let p = row.len();
        let mut xj = 1.;

        for j in 0..p {
            row.set(j, xj);
            xj *= t;
        }
    }

    fn solve_system(print_data: bool, t: MultilargeLinearType, c: &mut VecF64) {
        let mut w =
            MultilargeLinearWorkspace::new(t, P).expect("MultilargeLinearWorkspace::new failed");
        let mut x = MatF64::new(NROWS, P);
        let mut y = VecF64::new(NROWS);
        let mut r = Rng::new(RngType::default()).expect("Rng::new failed");

        let mut reg_param = VecF64::new(NLCURVE);
        let mut rho = VecF64::new(NLCURVE);
        let mut eta = VecF64::new(NLCURVE);

        let mut rowidx = 0;
        let mut t = 0.;

        while rowidx < N {
            // number of rows left to accumulate
            let nleft = N - rowidx;
            // number of rows in this block
            let nr = if NROWS > nleft { nleft } else { NROWS };

            let mut xv = x.submatrix(0, 0, nr, P);
            let mut yv = y.subvector(0, nr);

            // build (X,y) block with 'nr' rows
            for i in 0..nr {
                let mut row = xv.row(i);
                let fi = func(t);
                // noise
                let ei = r.gaussian(0.1 * fi);
                let yi = fi + ei;

                // construct this row of LS matrix
                build_row(t, &mut row);

                // set right hand side value with added noise
                yv.set(i, yi);

                if print_data && i % 100 == 0 {
                    println!("{} {}", t, yi);
                }

                t += DT;
            }

            // accumulate (X,y) block into LS system
            w.accumulate(&mut xv, &mut yv);

            rowidx += nr;
        }

        if print_data {
            println!();
            println!();
        }

        // compute L-curve
        if let Err(e) = w.lcurve(&mut reg_param, &mut rho, &mut eta) {
            eprintln!(
                "=== Method {} FAILED ===",
                w.name().expect("Failed to get name")
            );
            eprintln!("error: {:?}", e);
            return;
        }

        // solve large LS system and store solution in c
        let (rnorm, snorm) = w.solve(LAMBDA, c).unwrap();

        // compute reciprocal condition number
        let rcond = w.rcond().unwrap();

        eprintln!("=== Method {} ===\n", w.name().expect("Failed to get name"));
        eprintln!("condition number = {}", 1. / rcond);
        eprintln!("residual norm    = {}", rnorm);
        eprintln!("solution norm    = {}", snorm);

        // output L-curve
        for i in 0..NLCURVE {
            println!(
                "{:.12} {:.12} {:.12}",
                reg_param.get(i),
                rho.get(i),
                eta.get(i)
            );
        }
    }

    pub fn run() {
        let mut c_tsqr = VecF64::new(P);
        let mut c_normal = VecF64::new(P);

        // solve system with TSQR method
        solve_system(true, MultilargeLinearType::tsqr(), &mut c_tsqr);

        println!();
        println!();

        // solve system with Normal equations method
        solve_system(false, MultilargeLinearType::normal(), &mut c_normal);

        // output solution
        let mut v = VecF64::new(P);
        let mut t = 0.;

        while t <= 1. {
            let f_exact = func(t);
            build_row(t, &mut v);

            let f_tsqr = blas::d::dot(&v, &c_tsqr).unwrap();
            let f_normal = blas::d::dot(&v, &c_normal).unwrap();

            println!("{} {:.6} {:.6} {:.6}", t, f_exact, f_tsqr, f_normal);

            t += 0.01;
        }
    }
}

#[cfg(feature = "v2_2")]
fn main() {
    crate::example::run();
}

#[cfg(not(feature = "v2_2"))]
fn main() {
    eprintln!("You need to enable the `v2_2` feature to be able to run this example");
}
