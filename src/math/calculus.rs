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

fn fft(n: &usize, direction: Sign, input: &Vec<c64> ) -> Result<Vec<c64>, Box<dyn std::error::Error>> {

    let n = n.clone();
    let mut plan: C2CPlan64 = C2CPlan::aligned(&[n], direction, Flag::MEASURE).unwrap();
    let mut in_v = AlignedVec::new(n);
    let mut out_v= AlignedVec::new(n);

    for i in 0..n {
        in_v[i] = c64::new(input[i].re, input[i].im);
    }

    plan.c2c(&mut in_v, &mut out_v).unwrap();

    Ok(out_v.to_vec())
}

// TODO: NEEDS Refactoring. I suspect I was pressed for time

/// See thesis and references therein for an explanation of the split-step method.
/// Maybe should create a Hamiltonian for ffts as well
fn split_step_s3<'a>(system_size: &'a f64, deltat: &'a f64, omega: &'a f64, interaction_strength: &'a f64, data: &'a Vec<Complex64>) -> Result<Vec<Complex64>, Box<dyn std::error::Error>> {

    let q = PI * 40. / system_size.clone();
    let n = data.len();

    let mut k_vec: Vec<Complex64> = vec![Complex64::from(0f64); n]
        .iter()
        .enumerate()
        .map(|(idx, val)| -> Complex64 {
            if (idx as f64) < (n as f64)/2.0 {
                Complex64::from((TWO_PI / system_size) * (idx as f64))
            } else {
                Complex64::from(-1.0 *(TWO_PI / system_size) * ((n as f64) - (idx as f64)))
            }
        })
        .collect();

    let mut  k_fac: Vec<Complex64> = k_vec
        .iter()
        .map(|x| -> Complex64 {
            Complex64::exp(((-I * deltat * omega ) * x.pow(2) )/ 2f64 )
        })
        .collect();

    let mut ft = fft(&n, Sign::Forward, &data.to_vec()).unwrap();

    // Normalize
    let mult: Vec<Complex64> = ft.iter()
        .zip(k_fac.iter())
        .map(|(ft, k_fac )|{
            ft * k_fac
        })
        .collect();
    
    let normalized: Vec<Complex64> = mult
        .iter()
        .map(|el| { el / Complex64::from(n as f64)})
        .collect();

    let mut back_ft = fft(&n,Sign::Backward, &normalized).unwrap();

    let step = back_ft
        .iter()
        .enumerate()
        .map(|(idx, el)| {
            let xpos = (idx as f64) * system_size / (n as f64);
            let sinx = (q * xpos).sin();
            el * Complex64::exp(-I * omega * deltat * (interaction_strength * el.norm()) + potential(&xpos, &q, &true, &true) )
        })
        .collect();

    let transform_2 = fft(&n, Sign::Forward, &step).unwrap();

    // Normalize
    let mult_2: Vec<Complex64> = transform_2.iter()
        .zip(k_fac.iter())
        .map(|(ft, k_fac)| {
            ft * k_fac
        })
        .collect();

    let normalized_2: Vec<Complex64> = mult_2
        .iter()
        .map(|el| { el / Complex64::from(n as f64)})
        .collect();

    let out = fft(&n, Sign::Backward, &normalized_2).unwrap();

    Ok(out)
}

/// Not sure, but I think I never used this.
/// As stated in the thesis, It was just too slow.
fn split_step_s7<'a>(system_size: &'a f64, deltat: &'a f64, omega: &'a f64, interaction_strength: &'a f64, data: &'a Vec<Complex64>) -> Result<Vec<Complex64>, Box<dyn std::error::Error>>{

    //Bandrauk and Shen, J. Phys. A, 27:7147
    const OMEGA: [f64; 4] = [1.3151861047504085, -1.1776799841887, 0.2355733213359357, 0.784513610477560];
    let mut out = data.clone();
    out = split_step_s3(&system_size, &deltat, &OMEGA[3], &interaction_strength, &out)?;
    out = split_step_s3(&system_size, &deltat, &OMEGA[2], &interaction_strength, &out)?;
    out = split_step_s3(&system_size, &deltat, &OMEGA[1], &interaction_strength, &out)?;
    out = split_step_s3(&system_size, &deltat, &OMEGA[0], &interaction_strength, &out)?;
    out = split_step_s3(&system_size, &deltat, &OMEGA[1], &interaction_strength, &out)?;
    out = split_step_s3(&system_size, &deltat, &OMEGA[2], &interaction_strength, &out)?;
    out = split_step_s3(&system_size, &deltat, &OMEGA[3], &interaction_strength, &out)?;

    Ok(out)
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