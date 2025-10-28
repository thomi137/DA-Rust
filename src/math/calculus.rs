use fftw::*;
use fftw::array::AlignedVec;
use fftw::plan::{C2CPlan, C2CPlan64};
use fftw::types::{c64, Flag, Sign};
use num::{Complex, pow};
use num::complex::Complex64;
use num::pow::Pow;
use crate::physics::{PI, potential};

const TWO_PI: f64 = 2.0 * PI;
const I: Complex64 = Complex::I;

pub struct FftHelper {
    n: usize,
    plan_fwd: C2CPlan64,
    plan_bwd: C2CPlan64,
    buffer_in: AlignedVec<c64>,
    buffer_out: AlignedVec<c64>,
}

impl FftHelper {
    pub fn new(n: usize) -> Result<Self, Box<dyn std::error::Error>> {
        Ok(Self {
            n,
            plan_fwd: C2CPlan::aligned(&[n], Sign::Forward, Flag::MEASURE)?,
            plan_bwd: C2CPlan::aligned(&[n], Sign::Backward, Flag::MEASURE)?,
            buffer_in: AlignedVec::new(n),
            buffer_out: AlignedVec::new(n),
        })
    }

    pub fn fft(&mut self, direction: Sign, input: &[c64]) -> Result<Vec<c64>, Box<dyn std::error::Error>> {
        // Copy input data into FFTW-aligned buffer
        for (i, val) in input.iter().enumerate() {
            self.buffer_in[i] = *val;
        }

        // Run appropriate plan
        match direction {
            Sign::Forward => self.plan_fwd.c2c(&mut self.buffer_in, &mut self.buffer_out)?,
            Sign::Backward => self.plan_bwd.c2c(&mut self.buffer_in, &mut self.buffer_out)?,
        }

        Ok(self.buffer_out.to_vec())
    }
}


// TODO: NEEDS Refactoring. I suspect I was pressed for time

/// See thesis and references therein for an explanation of the split-step method.
/// Maybe should create a Hamiltonian for ffts as well
fn split_step_s3<'a>(
    fftw: &'a mut FftHelper,
    system_size: &'a f64,
    deltat: &'a f64,
    omega: &'a f64,
    interaction_strength: &'a f64,
    data: &[Complex64]) -> Result<Vec<Complex64>, Box<dyn std::error::Error>> {

    let q = PI * 40. / system_size.clone();
    let n = data.len();

    let mut k_vec: Vec<Complex64> = (0..n)
        .map(|idx| -> Complex64 {
            if (idx as f64) < (n as f64)/2.0 {
                Complex64::from((TWO_PI / system_size) * (idx as f64))
            } else {
                Complex64::from(-1.0 *(TWO_PI / system_size) * ((n as f64) - (idx as f64)))
            }
        })
        .collect();

    let mut k_fac: Vec<Complex64> = k_vec
        .iter()
        .map(|x| -> Complex64 {
            Complex64::exp(((-I * deltat * omega ) * x.pow(2) )/ 2f64 )
        })
        .collect();

    let mut ft = fftw.fft(Sign::Forward, data)?;

    for i in 0..n {
        ft[i] *= k_fac[i] / (n as f64);
    }

    let mut psi_x = fftw.fft(Sign::Backward, &ft)?;

    for (i, psi) in psi_x.iter_mut().enumerate() {
        let xpos = system_size * 0.5 - (i as f64) * system_size / (n as f64);
        let v = potential(&xpos, &q, &true, &true);
        let interaction = interaction_strength * psi.norm_sqr();
        let phase = -omega * deltat * (interaction + v);
        *psi *= Complex64::new(0.0, phase).exp();
    }

    let mut ft2 = fftw.fft(Sign::Forward, &psi_x)?;
    for i in 0..n {
        ft2[i] *= k_fac[i] / (n as f64);
    }
    let psi_out = fftw.fft(Sign::Backward, &ft2)?;


    let norm: f64 = psi_out.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
    let psi_out: Vec<c64> = psi_out.iter().map(|z| z / norm).collect();

    Ok(psi_out)
}

/// Not sure, but I think I never used this.
/// As stated in the thesis, It was just too slow.
fn split_step_s7<'a>(fftw:&'a mut FftHelper, system_size: &'a f64, deltat: &'a f64, omega: &'a f64, interaction_strength: &'a f64, data: &'a Vec<Complex64>) -> Result<Vec<Complex64>, Box<dyn std::error::Error>>{

    //Bandrauk and Shen, J. Phys. A, 27:7147
    const S7_COEFFS: [f64; 7] = [
        0.784513610477560,
        0.2355733213359357,
        -1.1776799841887,
        1.3151861047504085,
        -1.1776799841887,
        0.2355733213359357,
        0.784513610477560,
    ];

    let mut psi = data.to_vec();

    for &c in S7_COEFFS.iter(){
        let sub_dt = deltat * c;
        psi = split_step_s3(fftw, &system_size, &sub_dt, &omega, &interaction_strength, &psi)?;
    }

    Ok(psi)
}


#[cfg(test)]
mod tests {
    use fftw::types::c64;
    use crate::physics::FRAC_ROOT_TWO_PI;
    use super::*;

    #[test]
    fn fft_test() {
        let n: usize = 5;
        let input = vec![c64::new(FRAC_ROOT_TWO_PI, 0.0); n];

        let out = fft(&n, Sign::Forward, &input).unwrap();

        //  assert_eq!(out[0], c64::new(1f64, 0f64));
    }


}