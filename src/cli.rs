use clap::{arg, Parser};
use serde::{Serialize, Deserialize};

use crate::math::{EigenConfig, SplConfig, Jobz, Uplo};

#[derive(Parser, Debug)]
#[command(
author="Thomas Prosser <thomas@prosser.ch>",
version,
about = "Simulation of kicked Bose-Einstein Condensates",
long_about = None)]
pub struct Cli {

    #[arg(short, long)]
    pub config: Option<String>,

    #[command(flatten)]
    pub global: GlobalConfig,

    #[arg(long, value_parser = ["eig", "spl"])]
    pub alg: String,

    #[arg(long, value_enum, required_if_eq("alg", "eig"), default_value_t = Jobz::WithEigenvectors)]
    pub mode: Jobz,

    #[arg(long, value_enum, required_if_eq("alg", "eig"), default_value_t = Uplo::UpperTriangle)]
    pub symmetry: Uplo,

    #[arg(long, default_value_t = 0.1, required_if_eq("alg", "spl"))]
    pub dt: f64,

    #[arg(long, default_value_t = 1.0, required_if_eq("alg", "spl"))]
    pub omega: f64,
}

/// Learned from ChatGPT Seems Serde does not look at Clap
/// So default values decoupled from serializing/deserializing
fn default_step_num() -> usize { 1024 }
fn default_system_size() -> f64 { 10.0 }
fn default_interaction_strength() -> f64 { 0.0 }

#[derive(Parser, Debug, Clone, Serialize, Deserialize, Default)]
pub struct GlobalConfig {

    /// Optional output file name (required if flag is used)
    #[arg(short, long, requires="format", value_name = "FILE")]
    pub output: Option<String>,

    /// File format. Only mandatory if output is set
    #[arg(short, long, value_parser=["png", "svg"])]
    pub format: Option<String>,

    /// Number of steps to use
    #[arg(short='n', long, default_value_t=1024)]
    #[serde(default = "default_step_num")]
    pub step_num: usize,

    /// Length of the system
    #[arg(short, long, default_value_t=10.0)]
    #[serde(default = "default_system_size")]
    pub system_size: f64,

    /// Interaction strength of the bosons in the BEC, defaults to 0.0
    #[arg(short = 'g', long, default_value_t=0.0)]
    #[serde(default = "default_interaction_strength")]
    pub interaction_strength: f64,

    /// Puts the BEC into an optical trap, defaults to false
    #[arg(short, long, default_value_t=false)]
    pub trap: bool,

    /// Adds a harmonic optical lattice, defaults to false
    #[arg(short, long, default_value_t=false)]
    pub lattice: bool,

}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "alg", rename_all = "lowercase")]
pub enum AlgorithmConfig {
    Eig(EigenConfig),
    Spl(SplConfig),
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FullConfig {
    pub global: GlobalConfig,
    pub algorithm: AlgorithmConfig,
}

