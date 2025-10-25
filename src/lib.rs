pub mod linalg;
pub mod physics;
pub mod cli;

// use fftw::{ types::{ Sign, c64}, plan,};
// use num::{Complex};

// Type Definitions for pure convenience (and for making) public maybe.
// type Direction = Sign;
// type Plan = plan::C2CPlan64;

// Some useful Constants

//pub const E: f64 = std::f64::consts::E;

//const I: c64 = Complex::I;


#[cfg(test)]
mod tests {
    use crate::{
        linalg::{EigenConfig, Jobz, Uplo, solvers},
        physics::*,
    };

    #[test]
    fn test_eigenconfig() {
        let config = EigenConfig::init(
            Jobz::EigenValuesOnly,
            Uplo::UpperTriangle,
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

