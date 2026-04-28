//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{EigenNonSymmetricVWorkspace, MatC64, MatF64, VecC64, eigen};

fn main() -> Result<(), rgsl::Error> {
    let data = &mut [
        -1., 1., -1., 1., -8., 4., -2., 1., 27., 9., 3., 1., 64., 16., 4., 1.,
    ];
    let mut m = MatF64::from_mut_slice(data, 4, 4);
    let mut eval = VecC64::new(4);
    let mut evec = MatC64::new(4, 4);
    let mut w =
        EigenNonSymmetricVWorkspace::new(4).expect("EigenNonSymmetricVWorkspace::new failed...");

    w.nonsymmv(&mut m, &mut eval, &mut evec)?;

    eigen::nonsymmv_sort(&mut eval, &mut evec, eigen::Sort::AbsDesc)?;

    for i in 0..4 {
        let eval_i = eval.get(i);
        let evec_i = evec.column(i);
        println!("eigenvalue: {eval_i}");
        println!("eigenvector: ");
        for j in 0..4 {
            println!("  {}", evec_i.get(j));
        }
    }
    Ok(())
}
