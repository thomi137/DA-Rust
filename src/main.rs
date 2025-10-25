extern crate core;

use plotters::prelude::*;
use bec_rust::{
    linalg::{EigenConfig, Jobz, Uplo},
};
use bec_rust::physics::TridiagHamiltonian;
use bec_rust::solvers::tridiag_eigensolver;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    // these will be Command Line Arguments
    const NUM_STEPS: usize = 1024;
    let config = EigenConfig::init(
        Jobz::WithEigenvectors,
        Uplo::LowerTriangle,
        NUM_STEPS,
        10.
    );

    // let hamiltonian = Hamiltonian::new(&config,0.1,true, true);
    // let result = symmetric_eigensolver(&config, &hamiltonian.operator);

    let hamiltonian = TridiagHamiltonian::new(&config,0.1,true, true);
    let result = tridiag_eigensolver(&config, hamiltonian.vectors.diag, hamiltonian.vectors.offdiag );

    let eigenvectors = match result {
        Ok(result) => result.1,
        Err(msg) => panic!("{:#?}", msg)
    };

    let plotvec: Vec<_> = eigenvectors[0..NUM_STEPS]
        .iter()
        .enumerate()
        .map(|(idx, val)| {
            let xpos =  (10.0 * 0.5) - ((idx as f64) * 10.0)/(NUM_STEPS as f64);
            (xpos, val.abs())
        })
        .collect();


    let root = BitMapBackend::new("./plots/0.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("BEC Groundstate", ("sans-serif", 50).into_font())
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-5f64..5f64, 0f64..0.1f64)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            &mut plotvec.into_iter(),
            &BLUE,
        ))?;
       // .label("y = x^2")
       // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    //chart
    //    .configure_series_labels()
    //    .background_style(&WHITE.mix(0.8))
    //    .border_style(&BLACK)
    //    .draw()?;

    root.present()?;

    Ok(())
}