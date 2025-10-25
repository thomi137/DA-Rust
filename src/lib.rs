// use fftw::{ types::{ Sign, c64}, plan,};
// use num::{Complex};

// Type Definitions for pure convenience (and for making) public maybe.
// type Direction = Sign;
// type Plan = plan::C2CPlan64;

// Some useful Constants
pub const PI: f64 = std::f64::consts::PI;
pub const FRAC_ROOT_TWO_PI: f64 = 0.398942280401432677939946059934381868_f64;
//pub const E: f64 = std::f64::consts::E;

//const I: c64 = Complex::I;

pub mod physics {

    use num::pow;
    use crate::FRAC_ROOT_TWO_PI;
    use crate::linalg::EigenConfig;

    #[derive(Clone, Debug)]
    pub struct Hamiltonian {
        pub operator: Vec<f64>,
        pub interaction_strength: f64,
        pub trap: bool,
        pub lattice: bool,
    }
    impl Hamiltonian {
        pub fn new(config: &EigenConfig, interaction_strength: f64, trap: bool, lattice: bool) -> Hamiltonian {
            let config = config.clone();
            let num_steps = config.n as usize;
            let fnum_steps= num_steps as f64;
            let system_width = config.system_width;
            let operator = init_operator(&num_steps, &system_width, &fnum_steps, &interaction_strength, &trap, &lattice, |row, col| {
                let step_size = system_width / fnum_steps;
                if row == col {
                    let wave_number = ( crate::PI * 40. ) / system_width.clone();
                    let xpos = position(&col, &system_width, &fnum_steps);
                    let val = 1./(&step_size * &step_size)
                        + interaction_strength * pow(FRAC_ROOT_TWO_PI * f64::exp(- pow(xpos, 2) / 2.), 2)
                        + potential(&xpos, &wave_number, &lattice, &trap);
                    return val
                } else if (col as isize - row as isize).abs() == 1  {
                    return -0.5 * 1./(pow(step_size, 2))
                }
                0.0
            });

            Hamiltonian{ operator, interaction_strength, trap, lattice }
        }



    } // Impl Hamiltonian

    #[derive(Clone, Debug)]
    pub struct TridiagOperator {
        pub diag: Vec<f64>,
        pub offdiag: Vec<f64>,
    }

    #[derive(Clone, Debug)]
    pub struct TridiagHamiltonian {
        pub vectors: TridiagOperator,
        pub interaction_strength: f64,
        pub trap: bool,
        pub lattice: bool,
    }
    impl TridiagHamiltonian {
        pub fn new(config: &EigenConfig, interaction_strength: f64, trap: bool, lattice: bool) -> TridiagHamiltonian {
            let config = config.clone();
            let n = config.n as usize;
            let system_width = config.system_width;
            let fnum_steps= n as f64;
            let step_size = system_width / fnum_steps;

            let d = vec![0.0; n]
                .iter()
                .enumerate()
                .map(|(idx, _)| {
                    let wave_number = ( crate::PI * 40. ) / system_width.clone();
                    let xpos = position(&idx, &system_width, &fnum_steps);
                     1./(&step_size * &step_size)
                        + interaction_strength.clone() * pow(FRAC_ROOT_TWO_PI * f64::exp(- pow(xpos, 2) / 2.), 2)
                        + potential(&xpos, &wave_number, &lattice, &trap)

                })
                .collect();

            let e = vec![0.0; &n-1]
                .iter()
                .map(|_|{
                    -0.5 * 1./(pow(step_size, 2))
                }).collect();

            let vectors = TridiagOperator{diag: d, offdiag: e };

            TridiagHamiltonian{ vectors, interaction_strength, trap, lattice}
        }
    }


    pub fn init_operator(num_steps: &usize, _system_width: &f64, _fnum_steps: &f64, _interaction_strength: &f64, _lattice: &bool, _trap: &bool, f: impl Fn(usize, usize) -> f64) -> Vec<f64> {
        let mut operator: Vec<f64> = Vec::new();
        for row in 0..*num_steps {
            for col in 0..*num_steps {
                let row = row.clone();
                operator.push( f(row, col));
            }
        };
        operator
    }

    pub fn position(index: &usize, system_width: &f64, fnum_steps: &f64) -> f64 {
        let idx = (*index).clone();
        (system_width * 0.5) - (idx as f64) * (system_width / fnum_steps)
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
                let pot = 0.5 * location * location + 0.5 * wave_number * wave_number * &sinx_sq;
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
    impl Jobz{
        fn get_jobz(selection: Jobz) -> u8 {
            match selection {
                Jobz::EigenValuesOnly => b'N',
                Jobz::WithEigenvectors => b'V',
            }
        }
    }


    pub enum Uplo {
        UpperTriangle,
        LowerTriangle
    }
    impl Uplo {
        fn get_uplo(upper_or_lower: Uplo) -> u8 {
            match upper_or_lower {
                Uplo::UpperTriangle => b'U',
                Uplo::LowerTriangle => b'L',
            }
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
    #[derive(Clone)]
    pub struct EigenConfig {
        pub jobz: u8,
        pub uplo: u8,
        pub n: i32,
        pub system_width: f64,
        pub lda: i32,
        pub lwork: i32
    }
    impl EigenConfig {
        pub fn init(job: Jobz, upper_lower: Uplo, step_number: usize, system_width: f64 ) -> EigenConfig {
            let jobz = Jobz::get_jobz(job);
            let uplo = Uplo::get_uplo(upper_lower);
            let n = step_number as i32;
            let n_to_move = n.clone();
            let lda = n.clone();
            let lwork = 3*n-1;
            let system_width = system_width;

            EigenConfig { jobz, uplo, n: n_to_move, system_width, lwork, lda }
        }
    } // Impl Eigenconfig
} // mod Linalg

pub mod solvers {

    use lapack::*;
    use spinoff::{Color, Spinner, spinners};
    use crate::linalg::{EigenConfig};

    pub fn symmetric_eigensolver(config: &EigenConfig, hamiltonian: &Vec<f64>) -> Result<(Vec<f64>, Vec<f64>), String> {

        let cfg = config.clone();
        let rng = cfg.n.clone() as usize;
        let lwork = cfg.lwork.clone();
        let mut matrix = hamiltonian.clone();
        let mut w = vec![0.0; rng];
        let mut work = vec![0.0; cfg.lwork as usize];
        let mut info = 2939;

        unsafe {
            dsyev(cfg.jobz, cfg.uplo, cfg.n, &mut matrix, cfg.lda, &mut w, &mut work, lwork, &mut info);
        }

        if info == 0 {
           return if cfg.jobz ==  b'N' { Ok((w, vec![])) }
            else if cfg.jobz == b'V' { Ok((w, matrix))}
            else { Ok((vec![], vec![])) }
        }

        Err(format!("Algorithm failed: info is {info}, please consult https://www.netlib.org/lapack/explore-html-3.6.1/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html"))
    }

    pub fn tridiag_eigensolver(config: &EigenConfig, mut diag: Vec<f64>, mut offdiag: Vec<f64>) -> Result<(Vec<f64>, Vec<f64>), String> {
        let config = config.clone();
        let n = diag.len() as i32;
        let jobz = config.jobz;

        let compute_vectors = jobz == b'V';
        let mut z = if jobz == b'V' {
            vec![0.0f64; (n.clone() * n.clone()) as usize]
        } else {
            vec![0.0f64; 1]
        };

        let ldz = if compute_vectors { n.clone() } else { 1 };
        let mut work = vec![0.0f64; if compute_vectors { (2 * n.clone() - 2).max(1) as usize } else { 1 }];

        let mut info = 2934;

        let mut sp = Spinner::new(spinners::Aesthetic, "Calculating Eigenvalues/Eigenvectors...", Color::Cyan);
        unsafe {
            dstev(
                jobz,
                n,
                &mut diag,
                &mut offdiag,
                &mut z,
                ldz,
                &mut work,
                &mut info,
            );
        }


        if info == 0 {
            sp.success("Done calculating");
            return if !compute_vectors { Ok((diag, vec![])) }
            else if compute_vectors { Ok((diag, z))}
            else { Ok((vec![], vec![])) }
        }

        sp.fail("There has been an error");
        Err(format!("Algorithm failed: info is {}, please consult https://www.netlib.org/lapack/explore-html/db/daa/group__stev_ga774cb10a577ccdd97546ef1c5599f8ee.html#ga774cb10a577ccdd97546ef1c5599f8ee", info))
    }

} // mod Solvers

#[cfg(test)]
mod tests {
    use crate::linalg::EigenConfig;
    use crate::physics::potential;
    use super::{linalg, physics, solvers};

    #[test]
    fn test_eigenconfig() {
        let config = EigenConfig::init(
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
        let pos = physics::position(&5, &10., &10.);
        let pos2 = physics::position(&10, &10., &10.);

        assert_eq!(pos, 0.0);
        assert_eq!(pos2, -5.0);
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

    #[test]
    fn eigensolver_test(){
        let config = EigenConfig{ jobz: b'V', uplo: b'U', n: 3, lda: 3, lwork:8, system_width: 10.};
        let hamiltonian = vec![3.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 3.0];

        let res = solvers::symmetric_eigensolver(&config, &hamiltonian).unwrap();

        assert_eq!(res.0.len(), 3);
        assert_eq!(res.1.len(), 9);

        // from docs.
        for (one, another) in res.0.iter().zip(&[2.0, 2.0, 5.0]) {
            assert!((one - another).abs() < 1e-14);
        }
    }

    #[test]
    fn eigenvector_eigensolver_test(){
        let config = EigenConfig{ jobz: b'V', uplo: b'U', n: 3, lda: 3, lwork:8, system_width: 10.};
        let hamiltonian = vec![3.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 3.0];

        let res = solvers::symmetric_eigensolver(&config, &hamiltonian).unwrap();
        let vecs = res.1;

        assert_eq!(vecs.len(), 9);
    }

    #[test]
    fn only_eigenvalues_eigensolver_test(){
        let config = EigenConfig{ jobz: b'N', uplo: b'U', n: 3, lda: 3, lwork:8, system_width: 10.};
        let hamiltonian = vec![3.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 3.0];

        let res = solvers::symmetric_eigensolver(&config, &hamiltonian).unwrap();

        assert_eq!(res.0.len(), 3);
        assert_eq!(res.1.len(), 0);

        // from docs.
        for (one, another) in res.0.iter().zip(&[2.0, 2.0, 5.0]) {
            assert!((one - another).abs() < 1e-14);
        }
    }

    #[test]
    fn tridiag_eigensolver_test(){
        let config = EigenConfig{ jobz: b'V', uplo: b'U', n: 3, lda: 3, lwork:8, system_width: 10.};
        let diag: Vec<f64> = (1..=5).map(|x| { x as f64}).collect();
        let offdiag: Vec<f64> = vec![2.0; 4];
        let res = solvers::tridiag_eigensolver(&config, diag, offdiag).unwrap();

        assert_eq!(res.0.len(), 5);
        assert_eq!(res.1.len(), 25);

        let exp_ev = vec![ (-1.12) ,1.0 ,3.0, 5.0, 7.12 ];
        for (one, another) in res.0.iter().zip(&exp_ev) {
            assert!((one - another).abs() < 1e-2);
        }

    }

    #[test]
    fn tridiag_eigensolver_test_only_eigenvalues(){
        let config = EigenConfig{ jobz: b'N', uplo: b'U', n: 3, lda: 3, lwork:8, system_width: 10.};
        let diag: Vec<f64> = (1..=5).map(|x| { x as f64}).collect();
        let offdiag: Vec<f64> = vec![2.0; 4];
        let res = solvers::tridiag_eigensolver(&config, diag, offdiag).unwrap();

        assert_eq!(res.0.len(), 5);
        assert_eq!(res.1.len(), 0);

    }
}

