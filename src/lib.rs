/// Broken down into modules
pub mod math;
pub mod physics;
pub mod cli;

use std::fs;
use toml;

use cli::{GlobalConfig};
use crate::cli::{AlgorithmConfig, Cli, FullConfig};
use crate::math::{LapackConfig, SplConfig};

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

pub fn merge_globals(cli: &GlobalConfig, file: Option<GlobalConfig>) -> GlobalConfig {
    match file {
        Some(file_cfg) => GlobalConfig {
            output: cli.output.clone().or(file_cfg.output),
            format: cli.format.clone().or(file_cfg.format),
            step_num: if cli.step_num != 1024 { cli.step_num } else { file_cfg.step_num },
            system_size: if cli.system_size != 10.0 { cli.system_size } else { file_cfg.system_size },
            interaction_strength: if cli.interaction_strength != 0.0 {
                cli.interaction_strength
            } else {
                file_cfg.interaction_strength
            },
            trap: if cli.trap { true } else { file_cfg.trap },
            lattice: if cli.lattice { true } else { file_cfg.lattice },
        },
        None => cli.clone()
    }
}

pub fn merge_algorithm(cli: &Cli, file: Option<AlgorithmConfig>, globals: &GlobalConfig) -> AlgorithmConfig {
    match cli.alg.as_str() {
        "eig" => {
            // Take CLI values first
            let mode = cli.mode;
            let symmetry = cli.symmetry;
            // Build final config
            let final_config = LapackConfig { mode, symmetry }.build(globals);
            AlgorithmConfig::Eig(final_config)
        },
        "spl" => {
            let dt = cli.dt;
            let omega = cli.omega;
            AlgorithmConfig::Spl(SplConfig { dt, omega }.build(globals))
        },
        _ => {
            // Fallback to file config if present
            file.unwrap_or_else(|| panic!("Unknown algorithm: {}", cli.alg))
        }
    }
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
