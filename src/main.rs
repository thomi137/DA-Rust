use std::error::Error;
use num::pow;
use lapack::*;
use plotters::prelude::*;
use bec_rust::*;
use bec_rust::linalg::{EigenConfig, Jobz, Uplo};
use bec_rust::physics::Hamiltonian;
use bec_rust::


fn main() -> Result<(), Box<dyn std::error::Error>> {

    // these will be Command Line Arguments
    const NUM_STEPS: usize = 2048;

    let config = EigenConfig::init(
        Jobz::WithEigenvectors,
        Uplo::UpperTriangle,
        NUM_STEPS,
        10.
    );

    let hamiltonian = Hamiltonian::new(config,1.0,true, true);


    //println!("info is: {info}");
    //println!("{:?}", a);
    let mut plotvec: Vec<_> = a[0..NUM_STEPS]
        .iter()
        .map(|item| { item.clone()})
        .enumerate()
        .map(|(idx, val)| {
            let xpos = (&system_width * 0.5) - (idx as f64) * (&system_width/&fnum_steps);
            (xpos, val)
        })
        .collect();

    println!("{:#?}", plotvec);

    let root = BitMapBackend::new("./0.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("BEC Groundstate", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-5f64..5f64, 0f64..0.075f64)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            &mut plotvec.into_iter(),
            &BLUE,
        ))?;
       // .label("y = x^2")
       // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}