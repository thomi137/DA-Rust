use std::error::Error;
use lapack::*;
use num::complex::Complex64;
use crate::math::calculus::{FftHelper, split_step_s3, split_step_s7};
use super::*;

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
        return if cfg.jobz == b'N' { Ok((w, vec![])) } else if cfg.jobz == b'V' { Ok((w, matrix)) } else { Ok((vec![], vec![])) }
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
        return if !compute_vectors { Ok((diag, vec![])) } else if compute_vectors { Ok((diag, z)) } else { Ok((vec![], vec![])) }
    }

    Err(format!("Algorithm failed: info is {}, please consult https://www.netlib.org/lapack/explore-html/db/daa/group__stev_ga774cb10a577ccdd97546ef1c5599f8ee.html#ga774cb10a577ccdd97546ef1c5599f8ee", info))
}

pub enum SplitStepAlgorithm{
    S3,
    S7
}

pub fn split_step_solver(n: usize, system_size: f64, omega: f64, dt: f64,  g: f64, alg:  SplitStepAlgorithm, data: &[Complex64]) -> Result<Vec<Complex64>, Box<dyn Error>>{

    let out = match alg {
        SplitStepAlgorithm::S3 => {
            let mut fft = FftHelper::new(n)?;
            split_step_s3(&mut fft, &system_size, &dt, &omega, &g, data)?
        },
        SplitStepAlgorithm::S7 => {
            let mut fft = FftHelper::new(n)?;
            split_step_s7(&mut fft, &system_size, &dt, &omega, &g, data)?
        }
    };

    Ok(out)
}




#[cfg(test)]
mod tests {
    use crate::math::{EigenConfig, solvers};

    #[test]
    fn eigensolver_test() {
        let config = EigenConfig { jobz: b'V', uplo: b'U', n: 3, lda: 3, lwork: 8, system_width: 10. };
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
    fn eigenvector_eigensolver_test() {
        let config = EigenConfig { jobz: b'V', uplo: b'U', n: 3, lda: 3, lwork: 8, system_width: 10. };
        let hamiltonian = vec![3.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 3.0];

        let res = solvers::symmetric_eigensolver(&config, &hamiltonian).unwrap();
        let vecs = res.1;

        assert_eq!(vecs.len(), 9);
    }

    #[test]
    fn only_eigenvalues_eigensolver_test() {
        let config = EigenConfig { jobz: b'N', uplo: b'U', n: 3, lda: 3, lwork: 8, system_width: 10. };
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
    fn tridiag_eigensolver_test() {
        let epsilon = 1e-2f64;
        let config = EigenConfig { jobz: b'V', uplo: b'U', n: 3, lda: 3, lwork: 8, system_width: 10. };
        let diag: Vec<f64> = (1..=5).map(|x| { x as f64 }).collect();
        let offdiag: Vec<f64> = vec![2.0; 4];
        let res = solvers::tridiag_eigensolver(&config, diag, offdiag).unwrap();

        assert_eq!(res.0.len(), 5);
        assert_eq!(res.1.len(), 25);

        let exp_ev = vec![(-1.12f64), 1.0f64, 3.0f64, 5.0f64, 7.12464];
        for (one, another) in res.0.iter().zip(&exp_ev) {
            assert!((one - another).abs() < epsilon);
        }
    }

    #[test]
    fn tridiag_eigensolver_test_only_eigenvalues() {
        let config = EigenConfig { jobz: b'N', uplo: b'U', n: 3, lda: 3, lwork: 8, system_width: 10. };
        let diag: Vec<f64> = (1..=5).map(|x| { x as f64 }).collect();
        let offdiag: Vec<f64> = vec![2.0; 4];
        let res = solvers::tridiag_eigensolver(&config, diag, offdiag).unwrap();

        assert_eq!(res.0.len(), 5);
        assert_eq!(res.1.len(), 0);
    }
}
