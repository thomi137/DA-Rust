use plotters::prelude::*;
use plotters::coord::Shift;

use spinoff::{Spinner, spinners, Color};

use bec_rust::{build_algorithm_config, load_config_from_file, merge_algorithm, merge_globals, save_config_to_file};
use clap::Parser;
use bec_rust::cli::*;
use bec_rust::math::{build_solver, SolverResult};
use bec_rust::physics::TridiagHamiltonian;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let cli = Cli::parse();

    let mut file_globals: Option<GlobalConfig> = None;
    let mut file_alg: Option<AlgorithmConfig> = None;

    if let Some(ref path) = cli.config {
        let file_config = load_config_from_file(&path)?;
        file_globals = Some(file_config.global);
        file_alg = Some(file_config.algorithm);
    }
    let globals = GlobalConfig { ..cli.global.clone() };

    let run_config = FullConfig{
        global: merge_globals(&globals, file_globals),
        algorithm: merge_algorithm(&cli, file_alg, &globals)
    };

    let toml_str = toml::to_string_pretty(&run_config)?;
    println!("{}", toml_str);

    let solver = build_solver(run_config.algorithm).unwrap();

    let mut sp = Spinner::new(spinners::Aesthetic, "Starting calculation.", Color::Cyan);

    sp.update_text("Eigenvalue/Vector calculation");
    let hamiltonian = TridiagHamiltonian::new(&run_config.global);

    let result = solver.run(&globals,hamiltonian.vectors.diag, hamiltonian.vectors.offdiag, None)?;
    match result {
        SolverResult::SplitStep(psi) => {
            // psi is your Vec<Complex<f64>>
            println!("Wavefunction: {:?}", psi);
        }
        SolverResult::Eigen(eigenvals, eigenvecs) => {
            // eigenvals is your Vec<f64>
            println!("Eigenvalues: {:?}", eigenvals);
        }
    }
    /*
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

    println!("{}", eigenvectors.len());

    sp.success("Done 😃");

 */

/*

if let Some(output) = globals.output {

    let filename: String = output;

    // Just to make very sure.
    let format: String = globals.format.unwrap_or("png".to_string());

    let plotvec: Vec<(f64, f64)> = eigenvectors[0..globals.step_num]
        .iter()
        .enumerate()
        .map(|(idx, val)| {
            let xpos =  (globals.system_size * 0.5) - ((idx as f64) * globals.system_size)/(globals.step_num as f64);
            (xpos, val.abs())
        })
        .collect();

    sp.update_text("Starting plot");
    make_plot(&format, &filename, &plotvec).expect("TODO: panic message");
}



   // .label("y = x^2")
   // .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

//chart
//    .configure_series_labels()
//    .background_style(&WHITE.mix(0.8))
//    .border_style(&BLACK)
//    .draw()?;

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
        let _ = draw_chart(&root, &plotvec);
        root.present()?;
    },
    other =>  {
        eprintln!("Unsupported format: {}. Use png or svg.", other);
        std::process::exit(1);
    }
}
*/
    Ok(())
}

fn draw_chart<DB: DrawingBackend>(root: &DrawingArea<DB, Shift>, data: &Vec<(f64, f64)>) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {
    root.fill(&WHITE)?;

    let (x_min, x_max) = data.iter().copied().map(|(x, _)| x).fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| (lo.min(x), hi.max(x)));
    let (y_min, y_max) = data.iter().copied().map(|(_, y)| y).fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), y| (lo.min(y), hi.max(y)));

    let mut chart = ChartBuilder::on(&root)
        .caption("BEC Groundstate", ("sans-serif", 32).into_font())
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart.configure_mesh().draw()?;

    chart.draw_series(LineSeries::new(
        data.iter().cloned(),
        &BLUE, )).expect("TODO: panic message");
    Ok(())
}