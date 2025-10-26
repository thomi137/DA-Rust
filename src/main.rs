use plotters::prelude::*;
use plotters::coord::Shift;

use spinoff::{Spinner, spinners, Color};

use bec_rust::{
    linalg::{EigenConfig, Jobz, Uplo},
};
use clap::Parser;
use bec_rust::cli::*;
use bec_rust::physics::TridiagHamiltonian;
use bec_rust::linalg::solvers::tridiag_eigensolver;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let args = Cli::parse();

    let mut sp = Spinner::new(spinners::Aesthetic, "Starting calculation.", Color::Cyan);
    let config = EigenConfig::init(
        args.mode as Jobz,
        Uplo::LowerTriangle,
        args.step_num,
        args.system_size
    );

    //let hamiltonian = Hamiltonian::new(&config,0.1,true, true);
    // let result = symmetric_eigensolver(&config, &hamiltonian.operator);

    sp.update_text("Eigenvalue/Vector calculation");
    let hamiltonian = TridiagHamiltonian::new(&config,0.1,args.trap, args.lattice);
    let result = tridiag_eigensolver(&config, hamiltonian.vectors.diag, hamiltonian.vectors.offdiag );

    let eigenvectors = match result {
        Ok(result) => {
            sp.update_text("Groundstate found");
            result.1
        },
        Err(msg) => {
            sp.fail(&format!("Something went wrong: {}", msg));
            panic!("{:#?}", msg)
        }
    };

    let plotvec: Vec<(f64, f64)> = eigenvectors[0..args.step_num]
        .iter()
        .enumerate()
        .map(|(idx, val)| {
            let xpos =  (args.system_size * 0.5) - ((idx as f64) * args.system_size)/(args.step_num as f64);
            (xpos, val.abs())
        })
        .collect();

    sp.update_text("Starting plog");
    make_plot("png", "result", &plotvec).expect("TODO: panic message");


       // .label("y = x^2")
       // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    //chart
    //    .configure_series_labels()
    //    .background_style(&WHITE.mix(0.8))
    //    .border_style(&BLACK)
    //    .draw()?;

    sp.success("Done 😃");
    Ok(())
}

fn make_plot(format: &str, filename: &str, plotvec: &Vec<(f64, f64)>) -> Result<(), Box<dyn std::error::Error>> {

    let path = format!("./plots/{}.{}", filename, format);

    match format {
        "png" => {
            let root = BitMapBackend::new(&path,(800, 600) ).into_drawing_area();
            let _ = draw_chart(&root, &plotvec);
            root.present()?;
        },
        "svg" => {
            let root = SVGBackend::new(&path,(800, 600) ).into_drawing_area();
            let _ = draw_chart(&root, plotvec);
            root.present()?;
        },
        other =>  {
            eprintln!("Unsupported format: {}. Use png or svg.", other);
            std::process::exit(1);
        }
    }

    Ok(())
}

fn draw_chart<DB: DrawingBackend>(root: &DrawingArea<DB, Shift>, data: &Vec<(f64, f64)>) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {
    root.fill(&WHITE)?;

    let (x_min, x_max) = data.iter().copied().map(|(x, _)| x).fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| (lo.min(x), hi.max(x)));
    let (y_min, y_max) = data.iter().copied().map(|(_, y)| y).fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), y| (lo.min(y), hi.max(y)));

    let mut chart = ChartBuilder::on(&root)
        .caption("BEC Groundstate", ("sans-serif", 50).into_font())
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart.configure_mesh().draw()?;

    chart.draw_series(LineSeries::new(
        data.iter().cloned(),
        &BLUE, )).expect("TODO: panic message");

    Ok(())
}