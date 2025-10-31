/// Broken down into modules
pub mod math;
pub mod physics;
pub mod cli;

use std::error::Error;
use std::fs;
use toml;

use cli::{GlobalConfig};
use crate::cli::{AlgorithmConfig, AlgorithmSubcommand, Cli, FullConfig};
use crate::math::{LapackConfig, SplConfig, Uplo, Jobz};

pub type SolverFn<Cfg> = fn(&Cfg, Vec<f64>, Vec<f64>) -> Result<(Vec<f64>, Vec<f64>), String>;

/// Kinda simulates an interface like in C# or Java
pub trait ConfigBuilder {
    type Output;
    fn build(&self, globals: &GlobalConfig) -> Self::Output;
}

pub fn build_algorithm_config(cli: &Cli, globals: &GlobalConfig) -> AlgorithmConfig {
    match &cli.alg {
        AlgorithmSubcommand::Eig(args) => {
            let eig = LapackConfig{ mode: args.mode, symmetry: args.symmetry}.build(globals);
            AlgorithmConfig::Eig(eig)
        },
        AlgorithmSubcommand::Spl(args) => AlgorithmConfig::Spl(
                SplConfig{dt: args.dt, omega: args.omega, imag_time: args.imag_time}
            )
    }
}

pub fn load_config_from_file(path: &str) -> Result<FullConfig, String> {
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

pub fn merge_configs(cli: &Cli, file: Option<FullConfig>) -> FullConfig {

    let global = if let Some(f) = file.as_ref() {
        GlobalConfig {
            output: cli.global.output.clone().or(f.global.output.clone()),
            format: cli.global.format.clone().or(f.global.format.clone()),
            step_num: if cli.global.step_num != 1024 { cli.global.step_num} else { f.global.step_num },
            system_size: if cli.global.system_size != 10.0 { cli.global.system_size } else { f.global.system_size },
            interaction_strength: if cli.global.interaction_strength != 0.0 {
                cli.global.interaction_strength
            } else {
                f.global.interaction_strength
            },
            trap: if cli.global.trap { true } else { f.global.trap },
            lattice: if cli.global.lattice { true } else { f.global.lattice },
        }
    } else { cli.global.clone() };

    let byte_map: fn(AlgorithmConfig) -> (Jobz, Uplo) = |alg: AlgorithmConfig| {
      match alg {
          AlgorithmConfig::Eig(cfg) => {
              (
                  Jobz::map_jobz(cfg.jobz),
                  Uplo::map_uplo(cfg.uplo)
              )
          },
          _ => (Jobz::WithEigenvectors, Uplo::UpperTriangle)
      }
    };

    let algorithm = {
        match &cli.alg {
            AlgorithmSubcommand::Eig(cli_args) => {
                let (mode, symmetry) = if let Some(f) = file.as_ref().map(|f|{ f.algorithm.clone() }) {
                    let ( f_mode, f_symmetry ) = byte_map(f);
                    (
                         if matches!(cli_args.mode, Jobz::WithEigenvectors) { cli_args.mode.clone() } else { f_mode},
                         if matches!(cli_args.symmetry, Uplo::UpperTriangle) { cli_args.symmetry.clone() } else { f_symmetry }
                    )
                } else {
                        (cli_args.mode.clone(), cli_args.symmetry.clone())
                };

                let eig = LapackConfig{ mode, symmetry }.build(&global);
                AlgorithmConfig::Eig(eig)
            },

            AlgorithmSubcommand::Spl(cli_args) => {
                let (dt, omega, imag_time) = if let Some(AlgorithmConfig::Spl(f)) = file.as_ref().map(|f| f.algorithm.clone()) {
                    (
                        if cli_args.dt != 0.1 { cli_args.dt } else { f.dt },
                        if cli_args.omega != 1.0 { cli_args.omega } else { f.omega },
                        if cli_args.imag_time { true } else { f.imag_time },
                    )
                } else {
                    (cli_args.dt, cli_args.omega, cli_args.imag_time)
                };
                AlgorithmConfig::Spl(SplConfig{ dt, omega, imag_time })
            }
        }
    };

    FullConfig { global, algorithm }
}

pub mod serde_ascii {
    use serde::{self, Deserialize, Deserializer, Serializer};

    pub fn serialize<S>(value: &u8, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: Serializer,
    {
        let c = *value as char;
        serializer.serialize_str(&c.to_string())
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<u8, D::Error>
        where
            D: Deserializer<'de>,
    {
        let s: String = Deserialize::deserialize(deserializer)?;
        if s.len() != 1 {
            return Err(serde::de::Error::custom("Expected single ASCII character"));
        }
        Ok(s.as_bytes()[0])
    }
}
