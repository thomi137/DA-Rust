use std::error::Error;
use fftw::array::AlignedVec;
use fftw::plan::{C2CPlan, C2CPlan64};
use fftw::types::{c64, Flag, Sign};
use num::{Complex};
use num::complex::{Complex64};

use serde::{Deserialize, Serialize};


use crate::physics::{PI, potential};
use crate::{GlobalConfig};

const TWO_PI: f64 = 2.0 * PI;
const I: Complex64 = Complex::I;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplConfig {
    pub dt: f64,
    pub omega: f64,
    pub imag_time: bool
}

pub struct FftHelper {
    plan_fwd: C2CPlan64,
    plan_bwd: C2CPlan64,
    buffer_in: AlignedVec<c64>,
    buffer_out: AlignedVec<c64>,
}

impl FftHelper {
    pub fn new(n: usize) -> Result<Self, Box<dyn std::error::Error>> {
        Ok(Self {
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


pub fn split_step_s3(
    fftw: &mut FftHelper,
    global: &GlobalConfig,
    conf: &SplConfig,
    psi: &[Complex64]) -> Result<Vec<Complex64>, Box<dyn Error>>{

    let system_size = &global.system_size;
    let g = &global.interaction_strength;
    let deltat = &conf.dt;
    let omega = &conf.omega;
    let trap = global.trap.clone();
    let lattice = global.lattice.clone();
    let imag_time = conf.imag_time;
    let psi_out = split_step_gen(fftw, &system_size, &deltat, &omega, g, trap, lattice, imag_time, psi)?;
    Ok(psi_out)
}

pub fn split_step_s3_imag(fftw: &mut FftHelper,
                          global: &GlobalConfig,
                          conf: &SplConfig,
                          psi: &[Complex64],
) -> Result<Vec<Complex64>, Box<dyn Error>>{
    let system_size = &global.system_size;
    let g = &global.interaction_strength;
    let deltat = &conf.dt;
    let omega = &conf.omega;
    let trap = global.trap.clone();
    let lattice = global.lattice.clone();
    let imag_time = conf.imag_time;
    let psi_out = split_step_gen(fftw, system_size, deltat, omega, g, trap, lattice, imag_time, psi)?;
    Ok(psi_out)
}

/// Not sure, but I think I never used this.
/// As stated in the thesis, It was just too slow.
pub fn split_step_s7(fftw: &mut FftHelper,
                         global: &GlobalConfig,
                         conf: &SplConfig,
                         psi: &[Complex64] )-> Result<Vec<Complex64>, Box<dyn std::error::Error>>{

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

    let mut psi = psi.to_vec();

    for &_c in S7_COEFFS.iter(){
        psi = split_step_s3(fftw,
                            &global,
                            &conf,
                            &psi)?;
    }

    Ok(psi)
}

/// Not sure, but I think I never used this.
/// As stated in the thesis, It was just too slow.
pub fn split_step_s7_imag( fftw: &mut FftHelper,
                           global: &GlobalConfig,
                           conf: &SplConfig,
                           psi: &[Complex64])-> Result<Vec<Complex64>, Box<dyn std::error::Error>>{

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

    let mut psi = psi.to_vec();

    for &c in S7_COEFFS.iter(){
        psi = split_step_s3_imag(fftw,
                                   &global,
                                   &conf,
                                   &psi)?;
    }

    Ok(psi)
}

// TODO: NEEDS Refactoring. I suspect I was pressed for time

/// See thesis and references therein for an explanation of the split-step method.
/// Maybe should create a Hamiltonian for ffts as well
fn split_step_gen(
    fftw: & mut FftHelper,
    system_size: &f64,
    deltat: &f64,
    omega: &f64,
    interaction_strength: &f64,
    trap: bool,
    lattice: bool,
    imag_time: bool,
    data: &[Complex64],
) -> Result<Vec<Complex64>, Box<dyn std::error::Error>> {

    let q = PI * 40. / system_size.clone();
    let n = data.len();

    let k_vec: Vec<Complex64> = (0..n)
        .map(|idx| -> Complex64 {
            let dk = TWO_PI / system_size;
            let k = if idx < n / 2 {
                dk * (idx as f64)
            } else {
                dk * ((idx as f64) - (n as f64))
            };
            Complex64::new(k, 0.0)
        })
        .collect();

    let k_fac: Vec<Complex64> = k_vec
        .iter()
        .map(|k| -> Complex64 {
            if !imag_time {
                (-I * omega * deltat * k.powi(2)/ 2f64).exp()
            } else {
                (-omega * deltat * k.powi(2) / 2f64).exp()
            }
        })
        .collect();

    let mut ft = fftw.fft(Sign::Forward, data)?;

    for i in 0..n {
        ft[i] *= k_fac[i];
    }

    let mut psi_x = fftw.fft(Sign::Backward, &ft)?;
    for v in psi_x.iter_mut(){
        *v /= n as f64;
    }

    for (i, psi) in psi_x.iter_mut().enumerate() {
        let xpos = -0.5 * system_size + (i as f64) * system_size / (n as f64);
        let v = potential(&xpos, &q, trap, lattice);
        let interaction = interaction_strength * psi.norm_sqr();
        let phase = if !imag_time {
            Complex64::new(0.0,-omega * deltat * (interaction + v))
        } else {
            Complex64::new(-deltat * omega * (interaction + v), 0.0)
        };

        *psi *= phase.exp();
    }

    let mut ft2 = fftw.fft(Sign::Forward, &psi_x)?;
    for i in 0..n {
        ft2[i] *= k_fac[i];
    }
    let mut psi_out = fftw.fft(Sign::Backward, &ft2)?;
    for v in psi_out.iter_mut() {
        *v /= n as f64;
    }

    let dx = system_size / (n as f64);
    let norm: f64 = (psi_out.iter().map(|z| z.norm_sqr()).sum::<f64>() * dx).sqrt();
    let psi_out: Vec<c64> = psi_out.iter().map(|z| z / norm).collect();

    Ok(psi_out)
}


/// Integrate |psi|^2 over the spatial grid using composite Simpson's rule.
/// `system_size` = total domain length.
/// Assumes uniform spacing.
pub fn simpson_norm(psi: &[Complex64], system_size: f64) -> f64 {
    let n = psi.len();
    if n < 2 {
        return psi.get(0).map_or(0.0, |z| z.norm_sqr() * system_size);
    }

    let dx = system_size / (n as f64);

    // ensure even number of intervals for classic Simpson’s rule
    let m = if (n - 1) % 2 == 1 { n - 1 } else { n };

    let mut sum = psi[0].norm_sqr() + psi[m - 1].norm_sqr();

    // 4 * odd indices
    for i in (1..m - 1).step_by(2) {
        sum += 4.0 * psi[i].norm_sqr();
    }
    // 2 * even indices
    for i in (2..m - 1).step_by(2) {
        sum += 2.0 * psi[i].norm_sqr();
    }

    let mut integral = dx / 3.0 * sum;

    // if n-1 is odd (odd number of intervals), add trapezoid for last pair
    if m != n {
        integral += 0.5 * dx * (psi[m - 1].norm_sqr() + psi[m].norm_sqr());
    }

    integral
}


#[cfg(test)]
mod tests {
    use fftw::types::c64;
    use crate::physics::FRAC_ROOT_TWO_PI;
    use super::*;

    #[test]
    #[ignore]
    fn fft_test() {
        let n: usize = 5;
        let input = vec![c64::new(FRAC_ROOT_TWO_PI, 0.0); n];

        let out = FftHelper::fft(n, Sign::Forward, &input).unwrap();

        //  assert_eq!(out[0], c64::new(1f64, 0f64));
    }
}