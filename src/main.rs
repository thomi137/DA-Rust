use std::error::Error;
use num::pow;
use lapack::*;
use plotters::prelude::*;
use bec_rust::*;
use bec_rust::linalg::{EigenConfig, Jobz, Uplo};


fn main() -> Result<(), Box<dyn std::error::Error>> {

    // these will be Command Line Arguments
    const NUM_STEPS: usize = 2048;

    let config = EigenConfig::init(
        Jobz::WithEigenvectors,
        Uplo::UpperTriangle,
        NUM_STEPS,
        10.
    );

    let fnum_steps= config.n as f64;
    let system_width = config.system_width;
    let interaction_strength  = 5.;
    let wave_number = ( bec_rust::PI * 40. ) / &system_width;
    let step_size = &system_width / fnum_steps;

    let mut eigen: Vec<f64> = Vec::with_capacity(NUM_STEPS);
    let mut data: Vec<f64> = Vec::with_capacity(NUM_STEPS);

    for idx in 0..NUM_STEPS {
        let xpos = (&system_width * 0.5) - (idx as f64 * &system_width)/&fnum_steps;
        data.push( FRAC_ROOT_TWO_PI * f64::exp(- pow(xpos, 2) / 2.) )
    }

    let mut a: Vec<f64> = Vec::new();
    for i in 0..NUM_STEPS{
        let xpos = (&system_width * 0.5) - (i as f64 * &system_width)/&fnum_steps;
        let sinx = f64::sin(&wave_number * xpos);
        let datsq = data[i] * data[i];
        for j in 0..NUM_STEPS{
            let val: f64 = if i == j {
                1./(&step_size * &step_size)
                    + &interaction_strength * &datsq
                    + 0.5 * (&xpos * &xpos)
                    + 0.5 * (&wave_number * & wave_number) * (&sinx * &sinx)
            } else if  i == j + 1  {
                    -0.5 * 1./(&step_size * &step_size)
            } else { 0.};
            a.push(val);
        }
    }
    let length = a.len();
    // println!("a has length: {length}, 0th element: {}", a[0]);

    let n = NUM_STEPS as i32;
    let lda = n;
    let lwork = 3*&n-1;
    let mut w = vec![0.0; n as usize];
    let mut work = vec![0.0; lwork as usize];
    let mut info = 2939;

    unsafe {
        dsyev(config.jobz, config.uplo, config.n, &mut a, lda, &mut w, &mut work, lwork, &mut info);
    }

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