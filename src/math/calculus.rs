use fftw::*;
use fftw::plan::{C2CPlan64};
use fftw::types::c64;
use num::Complex;

type Direction = types::Sign;
type Plan = C2CPlan64;

fn fft(&n: usize, direction: Direction, plan: Plan, &input: Vec<c64> ) -> Result<Vec<c64>, error::Error> {


    Ok(())
}