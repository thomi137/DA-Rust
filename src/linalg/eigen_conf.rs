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

#[cfg(test)]
mod tests {
    use crate::linalg::{EigenConfig, Jobz, Uplo};

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
}
