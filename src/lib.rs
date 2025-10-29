/// Broken down into modules
pub mod math;
pub mod physics;
pub mod cli;

use std::fs;
use toml;

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

fn load_config_from_file(path: &str) -> Result<GlobalConfig, String> {
    let toml_str = fs::read_to_string(path)
        .map_err(|e| format!("Error reading {}: {}", path, e))?;
    toml::from_str(&toml_str)
        .map_err(|e| format!("Error parsing {}: {}", path, e))
}

fn save_config_to_file(config: &GlobalConfig, path: &str) -> Result<(), String> {
    let toml_str = toml::to_string_pretty(config)
        .map_err(|e| format!("Error serializing: {}", e))?;
    fs::write(path, toml_str)
        .map_err(|e| format!("Error writing {}: {}", path, e))
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

/*
pub fn parse_config_file(path: &str) -> SimulationConfig {
    let content = fs::read_to_string(path).expect("Failed to read config file");
    serde_json::from_str(&content).expect("Failed to parse config file")
}

 */


