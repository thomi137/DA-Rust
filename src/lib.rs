/// Broken down into modules
pub mod math;
pub mod physics;
pub mod cli;

use std::fs;
use cli::{GlobalConfig};
use crate::cli::{AlgorithmConfig, Cli};
use crate::math::LapackConfig;

pub type SolverFn<Cfg> = fn(&Cfg, Vec<f64>, Vec<f64>) -> Result<(Vec<f64>, Vec<f64>), String>;

/// Kinda simulates an interface like in C# or Java
pub trait ConfigBuilder {
    type Output;
    fn build(&self, globals: &GlobalConfig) -> Self::Output;
}

pub fn build_algorithm_config(cli: &Cli, globals: &GlobalConfig) -> AlgorithmConfig {
    match cli.alg.as_str() {
        "eig" => {
            let eig = LapackConfig{ mode: cli.mode, symmetry: cli.symmetry}.build(globals);
            AlgorithmConfig::Eig(eig)
        },
        /*
        "spl" => {
            let spl = SplConfig::init(cli.dt, cli.omega, globals);
            AlgorithmConfig::Spl(spl)
        }*/
        _ => panic!("Unknown algorithm: {}", cli.alg),
    }
}

/*
pub fn parse_config_file(path: &str) -> SimulationConfig {
    let content = fs::read_to_string(path).expect("Failed to read config file");
    serde_json::from_str(&content).expect("Failed to parse config file")
}

 */


