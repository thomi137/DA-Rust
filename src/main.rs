use plotters::prelude::*;
use bec_rust::{
    linalg::{EigenConfig, Jobz, Uplo},
};
use clap::Parser;
use bec_rust::cli::*;
use bec_rust::physics::TridiagHamiltonian;
use bec_rust::linalg::solvers::tridiag_eigensolver;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let args = Cli::parse();

    let config = EigenConfig::init(
        args.mode as Jobz,
        Uplo::LowerTriangle,
        args.step_num,
        args.system_size
    );

    //let hamiltonian = Hamiltonian::new(&config,0.1,true, true);
    // let result = symmetric_eigensolver(&config, &hamiltonian.operator);

    let hamiltonian = TridiagHamiltonian::new(&config,0.1,args.trap, args.lattice);
    let result = tridiag_eigensolver(&config, hamiltonian.vectors.diag, hamiltonian.vectors.offdiag );

    let eigenvectors = match result {
        Ok(result) => result.1,
        Err(msg) => panic!("{:#?}", msg)
    };

    let plotvec: Vec<_> = eigenvectors[0..args.step_num]
        .iter()
        .enumerate()
        .map(|(idx, val)| {
            let xpos =  (args.system_size * 0.5) - ((idx as f64) * 10.0)/(args.step_num as f64);
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