use clap::*;
use serde::{Serialize, Deserialize};

use crate::{ConfigBuilder, GlobalConfig};

/// Linear Algebra used.
/// at the moment, the eigenvalues and eigenvectors of a tridiagonal Matrix are used,
/// so we will use LAPACK's dsyev
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Serialize, Deserialize)]
pub enum Jobz {
    #[clap(name="values")]
    #[serde(rename = "values")]
    EigenValuesOnly,

    #[clap(name="vectors")]
    #[serde(rename = "vectors")]
    WithEigenvectors
}
impl Jobz{
    fn get_jobz(&self) -> u8 {
        match self {
            Jobz::EigenValuesOnly => b'N',
            Jobz::WithEigenvectors => b'V',
        }
    }

    pub fn map_jobz(byte: u8) -> Jobz{
        match byte {
            b'N' => Jobz::EigenValuesOnly,
            b'V' => Jobz::WithEigenvectors,
            _ => Jobz::WithEigenvectors
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Serialize, Deserialize)]
pub enum Uplo {
    #[clap(name="upper")]
    #[serde(rename = "upper")]
    UpperTriangle,

    #[clap(name="lower")]
    #[serde(rename = "lower")]
    LowerTriangle
}
impl Uplo {
    fn get_uplo(&self) -> u8 {
        match self {
            Uplo::UpperTriangle => b'U',
            Uplo::LowerTriangle => b'L',
        }
    }

    pub fn map_uplo(byte: u8) -> Uplo {
        match byte {
            b'U' => Uplo::UpperTriangle,
            b'L' => Uplo::LowerTriangle,
            _ => Uplo::UpperTriangle
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LapackConfig {
    pub(crate) mode: Jobz,
    pub(crate) symmetry: Uplo
}
impl ConfigBuilder for LapackConfig {
    type Output = EigenConfig;

    fn build(&self, globals: &GlobalConfig) -> Self::Output {
        let jobz = self.mode.get_jobz();
        let uplo = self.symmetry.get_uplo();
        let lda = globals.clone().step_num as i32;
        let lwork = 3 * globals.step_num as i32 -1;

        EigenConfig{
            jobz,
            uplo,
            lda,
            lwork
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EigenConfig {
    pub jobz: u8,
    pub uplo: u8,
    pub lda: i32,
    pub lwork: i32
}

#[cfg(test)]
mod tests {
    use crate::math::{EigenConfig, Jobz, Uplo};

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
