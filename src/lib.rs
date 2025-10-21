use fftw::{
    types::{
        Sign,
        c64
    },
    plan,
};

use num::{Complex};
use lapack::*;

// Type Definitions for pure convenience (and for making) public maybe.
type Direction = Sign;
type Plan = plan::C2CPlan64;

// Somee useful Constants
pub const PI: f64 = std::f64::consts::PI;
pub const FRAC_ROOT_TWO_PI: f64 = 0.398942280401432677939946059934381868_f64;
pub const E: f64 = std::f64::consts::E;

const I: c64 = Complex::I;

pub mod physics {

    use num::pow;
    use crate::linalg::EigenConfig;

    pub struct Hamiltonian {
        pub pot: Vec<f64>,
        pub interaction_strength: f64,
    }
    impl Hamiltonian {
        pub fn new(config: EigenConfig, interaction_strength: f64, lattice: bool, trap: bool) -> Hamiltonian {
            let config = config;
            let fnum_steps= config.n as f64;
            let system_width = config.system_width;
            let wave_number = ( crate::PI * 40. ) / &system_width;
            let step_size = &system_width / fnum_steps;
            let mut pot: Vec<f64> = Vec::new();
            for idx in 0..config.n {
                let xpos = (&system_width * 0.5) - (idx as f64 * &system_width)/&fnum_steps;
                pot.push(potential(xpos, &wave_number, &lattice, &trap) );
            };
            Hamiltonian{ pot, interaction_strength }
        }
    }

    fn potential(location: f64, wave_number: &f64, trap: &bool, lattice: &bool) -> f64 {
        let sinx = f64::sin( wave_number * location );

        match (trap, lattice) {
            (true, false) => &location * &location + 0.5,
            (false, true) => {
                let pot = 0.5 * pow(sinx, 2) * &location * &location;
                pot
            },
            (true, true) => {
                let sinx_sq = pow(sinx, 2);
                let lattice = 0.5 * &sinx_sq * &location * &location;
                let pot = 0.5 * &sinx_sq * &location * &location + 0.5 * &sinx_sq * &location * &location;
                pot
            },
            _ => 0.
        }
    }
}

/// Linear Algebra used.
/// at the moment, the eigenvalues and eigenvectors of a tridiagonal Matrix are used,
/// so we will use LAPACK's dsyev
pub mod linalg {

    pub enum Jobz {
        EigenValuesOnly,
        WithEigenvectors
    }
    fn get_jobz(selection: Jobz) -> u8 {
        match selection {
            Jobz::EigenValuesOnly => b'N',
            Jobz::WithEigenvectors => b'V',
        }
    }

    pub enum Uplo {
        UpperTriangle,
        LowerTriangle
    }
    fn get_uplo(upper_or_lower: Uplo) -> u8 {
        match upper_or_lower {
            Uplo::UpperTriangle => b'U',
            Uplo::LowerTriangle => b'L',
        }
    }

    /// #Config
    ///
    /// Configuration struct for LAPACK Functions in General
    /// Parameters:
    /// * `job`: Eigenvalues only, or also eigenvectors.
    /// * `upper_lower`: Upper or lower part of symmetric matrix not 0.
    /// * `n`: Matrix order. usize. We only need one rank, since we want a symmetric matrix.
    /// * `system_size`: width of the system
    pub struct EigenConfig {
        pub jobz: u8,
        pub uplo: u8,
        pub n: i32,
        pub system_width: f64,
    }
    impl EigenConfig {
        pub fn init(job: Jobz, upper_lower: Uplo, step_number: usize, system_width: f64 ) -> EigenConfig {
            let jobz = get_jobz(job);
            let uplo = get_uplo(upper_lower);
            let n = step_number as i32;
            let system_width = system_width;

            EigenConfig { jobz, uplo, n, system_width }
        }
    }
}
