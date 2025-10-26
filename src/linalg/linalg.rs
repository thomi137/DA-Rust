use clap::ValueEnum;

/// Linear Algebra used.
/// at the moment, the eigenvalues and eigenvectors of a tridiagonal Matrix are used,
/// so we will use LAPACK's dsyev
    #[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
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

    #[derive(Clone, Debug)]
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
    #[derive(Clone, Debug)]
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

    pub mod solvers {
        use lapack::*;
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
    }
