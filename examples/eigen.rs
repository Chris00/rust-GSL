//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{EigenSymmetricVWorkspace, MatF64, VecF64, eigen};

fn main() -> Result<(), rgsl::Error> {
    let data = &mut [
        1.,
        1. / 2.,
        1. / 3.,
        1. / 4.,
        1. / 2.,
        1. / 3.,
        1. / 4.,
        1. / 5.,
        1. / 3.,
        1. / 4.,
        1. / 5.,
        1. / 6.,
        1. / 4.,
        1. / 5.,
        1. / 6.,
        1. / 7.,
    ];
    let mut m = MatF64::from_mut_slice(data, 4, 4);
    let mut eval = VecF64::new(4);
    let mut evec = MatF64::new(4, 4);
    let mut w = EigenSymmetricVWorkspace::new(4).expect("EigenSymmetricVWorkspace::new failed...");

    w.symmv(&mut m, &mut eval, &mut evec)?;

    eigen::symmv_sort(&mut eval, &mut evec, eigen::Sort::AbsAsc)?;

    for i in 0..4 {
        let eval_i = eval.get(i);
        let evec_i = evec.column(i);
        println!("eigenvalue: {}", eval_i);
        println!("eigenvector: {:?}", evec_i);
    }

    Ok(())
}
