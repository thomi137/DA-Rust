use fftw::*;
use fftw::plan::{C2CPlan64};
use fftw::types::c64;
use num::Complex;

type Direction = types::Sign;
type Plan = C2CPlan64;

fn fft(n: &usize, direction: &Direction, plan: Plan, input: &Vec<c64> ) -> Result<(), error::Error> {

    let input: Vec<f64>;
    let output: Vec<f64>;
    //let plan = C2CPlan64(&input, &output, Direction::Forward);;


    Ok(())
}