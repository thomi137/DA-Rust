use fftw::{
    types::{
        Sign,
        c64
    },
    plan,
};

use num::{Complex};

// Type Definitions for pure convenience (and for making) public maybe.
type Direction = Sign;
type Plan = plan::C2CPlan64;

// Some useful Constants
pub const PI: f64 = std::f64::consts::PI;
pub const FRAC_ROOT_TWO_PI: f64 = 0.398942280401432677939946059934381868_f64;
pub const E: f64 = std::f64::consts::E;

const I: c64 = Complex::I;

pub mod physics {

    use num::pow;
    use crate::FRAC_ROOT_TWO_PI;
    use crate::linalg::EigenConfig;

    pub struct Hamiltonian {
        pub operator: Vec<f64>,
        pub interaction_strength: f64,
        pub trap: bool,
        pub lattice: bool,
    }
    impl Hamiltonian {
        pub fn new(&mut self, config: EigenConfig, interaction_strength: f64, trap: bool, lattice: bool) -> Hamiltonian {
            let fnum_steps= config.n as f64;
            let system_width = config.system_width;
            let operator: Vec<f64> = Vec::new();

            self.init_operator(&config.n, &system_width, &fnum_steps, &interaction_strength, &trap, &lattice );

            Hamiltonian{ operator, interaction_strength, trap, lattice }
        }

        fn init_operator(&mut self, mut len: &i32, system_width: &f64, fnum_steps: &f64, interaction_strength: &f64, lattice: &bool, trap: &bool) {
            for row in 0..*len {
                for col in 0..*len {
                    if row == col {
                        let column = col.clone();
                        let step_size = system_width / fnum_steps;
                        let wave_number = ( crate::PI * 40. ) / system_width;
                        let xpos = position(&column, &system_width, &fnum_steps);
                        let val = 1./(&step_size * &step_size)
                            + interaction_strength * (FRAC_ROOT_TWO_PI * f64::exp(- pow(xpos, 2) / 2.))
                            + potential(&xpos, &wave_number, lattice, trap);
                        self.operator.push(val);
                    } else if  row == col + 1  {
                        let step_size = system_width / fnum_steps;
                        self.operator.push(-0.5 * 1./(pow(step_size, 2)));
                    } else { self.operator.push(0.) };
                }
            };
        }



    } // Impl Hamiltonian

    pub fn position(index: &i32, system_width: &f64, fnum_steps: &f64) -> f64 {
        let idx = (*index).clone();
        -(system_width * 0.5) + (idx as f64) * (system_width / fnum_steps)
    }

    pub fn potential(location: &f64, wave_number: &f64, trap: &bool, lattice: &bool) -> f64 {
        let sinx = f64::sin( wave_number * location );

        match (trap, lattice) {
            (true, false) => location * location * 0.5,
            (false, true) => {
                let pot = 0.5 * pow(sinx, 2) * location * location;
                pot
            },
            (true, true) => {
                let sinx_sq = pow(sinx, 2);
                let pot = 0.5 * &sinx_sq * location * location + 0.5 * &sinx_sq * location * location;
                pot
            },
            _ => 0.
        }
    }
} // Mod Physics

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
    } // Impl Eigenconfig
} // mod Linalg

#[cfg(test)]
mod tests {
    use crate::physics::{Hamiltonian, potential};
    use super::{linalg, physics};

    #[test]
    fn test_eigenconfig() {
        let config = linalg::EigenConfig::init(
            linalg::Jobz::EigenValuesOnly,
            linalg::Uplo::UpperTriangle,
            10,
            10.,
        );

        assert_eq!(config.jobz, b'N');
        assert_eq!(config.uplo, b'U');
    }

    #[test]
    fn test_position(){
        let pos = super::physics::position(&5, &10., &10.);
        let pos2 = physics::position(&10, &10., &10.);

        assert_eq!(pos, 0.0);
        assert_eq!(pos2, 5.0);
    }

    // Trap only, Harmonic Oscillator
    #[test]
    fn test_potential(){
        let osc_pot = potential(&0., &1., &true, &false );
        let osc_pot2 = potential(&-5.0, &1., &true, &false );

        assert_eq!(osc_pot, 0.);
        assert_eq!(osc_pot2, 12.5);
    }

    // Lattice only. Sin^2 -like potential
    #[test]
    fn test_lat_potential(){
        let lat_pot = potential(&0., &1., &false, &true );

        assert_eq!(lat_pot, 0.);
    }

}

