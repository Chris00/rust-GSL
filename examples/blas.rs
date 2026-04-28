//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{
    MatF64,
    blas::{self, Transpose},
};

fn main() -> Result<(), rgsl::Error> {
    let a0 = &mut [0.11, 0.12, 0.13, 0.21, 0.22, 0.23];
    let b0 = &mut [1011., 1012., 1021., 1022., 1031., 1032.];
    let c0 = &mut [0., 0., 0., 0.];

    let a = MatF64::from_mut_slice(a0, 2, 3);
    let b = MatF64::from_mut_slice(b0, 3, 2);
    let mut c = MatF64::from_mut_slice(c0, 2, 2);

    blas::d::gemm(
        Transpose::NoTranspose,
        Transpose::NoTranspose,
        1.,
        &a,
        &b,
        0.,
        &mut c,
    )?;

    drop(c);
    println!("[ {}, {}", c0[0], c0[1]);
    println!("  {}, {} ]", c0[2], c0[3]);
    Ok(())
}
