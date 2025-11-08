/// Broken down into modules
pub mod math;
pub mod physics;
pub mod cli;

use std::fs;

use serde::{Serialize, Deserialize};
use toml;

use crate::cli::{GlobalConfig, AlgorithmConfig, AlgorithmSubcommand, Cli, FullConfig};
use crate::math::{LapackConfig, SplConfig, Uplo, Jobz, EigenConfig};

pub type SolverFn<Cfg> = fn(&Cfg, Vec<f64>, Vec<f64>) -> Result<(Vec<f64>, Vec<f64>), String>;

/// Use default values and keep in sync with CLI
/// Global Values:
const DEFAULT_STEP_NUM: usize = 1024usize;
const DEFAULT_SYSTEM_SIZE: f64 = 10f64;
const DEFAULT_INTERACTION_STRENGTH: f64 = 0f64;

/// Eigensolver Default values
/// Note that UPLO is really only needed by LAPACK
/// if we are talking Symmetric Matrices. In a time-independent
/// System, we would expect the Hamiltonian to be represented by
/// a symmetric matrix because the observables need to be real numbers.
/// and the set of eigenvectors forms an orthogonal basis as a consequence of the
/// [Spectral Theorem](https://en.wikipedia.org/wiki/Spectral_theorem)
/// In our special case, we actually only need a tridiagonal Matrix since as per the Thesis,
/// the interaction is small an concerns nearest neighbours, so we can further speed up calcualtion
/// and memory by using DSTEV instead of DSYEV from LAPACK
const DEFAULT_JOBZ: Jobz = Jobz::WithEigenvectors;
const DEFAULT_UPLO: Uplo = Uplo::UpperTriangle;

/// Spl Default Values
/// Note that we do not need a default for a boolean imag_time,
/// which forces imaginary time propagation in the split-step by making
/// the Gaussian exponent real (and therefore not unitary)
///
/// A boolean flag is there or it is not there, so if somebody says it is true,
/// it is true, else if everyone keeps quiet, it is false.
const DEFAULT_DELTAT: f64  = 0.01f64;
const DEFAULT_OMEGA: f64 = 1f64;

/// Kinda simulates an interface like in C# or Java
pub trait ConfigBuilder {
    type Output;
    fn build(&self, globals: &GlobalConfig) -> Self::Output;
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FileConfig {
    pub global: Option<GlobalConfig>,
    pub eig: Option<EigenConfig>,
    pub spl: Option<SplConfig>,
}

pub fn load_config_from_file(path: &str) -> Result<FileConfig, String> {
    let toml_str = fs::read_to_string(path)
        .map_err(|e| format!("Error reading {}: {}", path, e))?;
    toml::from_str(&toml_str)
        .map_err(|e| format!("Error parsing {}: {}", path, e))
}


pub fn save_config_to_file(config: &FullConfig, path: &str) -> Result<(), String> {
    let toml_str = toml::to_string_pretty(&config)
        .map_err(|e| format!("Error serializing: {}", e))?;
    fs::write(path, toml_str)
        .map_err(|e| format!("Error writing {}: {}", path, e))
}


pub fn merge_configs(cli: &Cli, f_conf: Option<FileConfig>) -> FullConfig {

    let global = match f_conf.as_ref() {
        Some(file) => {
            let globals = file.global.clone().unwrap_or(cli.global.clone());
            get_global_config(&cli, globals)
        }
        None => cli.global.clone()
    };

    let algorithm = get_algorithm_config(&cli, f_conf, &global);

    FullConfig { global, algorithm }
}

fn get_global_config(cli: &Cli, file_g: GlobalConfig) -> GlobalConfig {

    GlobalConfig {
        output: cli.global.output.clone().or(file_g.output.clone()),
        format: cli.global.format.clone().or(file_g.format.clone()),
        step_num: if cli.global.step_num != DEFAULT_STEP_NUM { cli.global.step_num } else { file_g.step_num },
        system_size: if cli.global.system_size != DEFAULT_SYSTEM_SIZE { cli.global.system_size } else { file_g.system_size },
        interaction_strength: if cli.global.interaction_strength != DEFAULT_INTERACTION_STRENGTH {
            cli.global.interaction_strength
        } else {
            file_g.interaction_strength
        },
        trap: if cli.global.trap { true } else { file_g.trap },
        lattice: if cli.global.lattice { true } else { file_g.lattice },
    }

}

/// Gets config for different algorithms. Needs to take into account
/// that certain algorithms derive some parameters according to their input
fn get_algorithm_config(cli: &Cli, f_conf: Option<FileConfig>, global: &GlobalConfig) -> AlgorithmConfig {
    match &cli.alg {
        AlgorithmSubcommand::Eig(cli_args) => {
            let e_conf = f_conf.and_then(|file| file.eig);
            let (mode, symmetry) = e_conf
                .and_then(|v| Some((Jobz::map_jobz(v.jobz), Uplo::map_uplo(v.uplo))))
                .map(|(j, u)| {
                    (
                        if matches!(cli_args.mode, Jobz::WithEigenvectors) { cli_args.mode.clone() } else { j },
                        if matches!(cli_args.symmetry, Uplo::UpperTriangle) { cli_args.symmetry.clone() } else { u }
                    )
                })
                .unwrap_or((Jobz::WithEigenvectors, Uplo::UpperTriangle));

            let eig = LapackConfig { mode, symmetry }.build(&global);
            AlgorithmConfig::Eig(eig)
        },

        AlgorithmSubcommand::Spl(cli_args) => {
            let s_conf = f_conf.and_then(|file| file.spl);
            let (dt, omega, imag_time) = s_conf.map(|s| {
                (
                    if cli_args.dt != DEFAULT_DELTAT { cli_args.dt.clone() } else { s.dt },
                    if cli_args.omega != DEFAULT_OMEGA { cli_args.omega.clone() } else { s.omega },
                    if cli_args.imag_time { true } else { s.imag_time },
                )
            }).unwrap_or((cli_args.dt, cli_args.omega, cli_args.imag_time));

            AlgorithmConfig::Spl(SplConfig { dt, omega, imag_time })
        }
    }
}
