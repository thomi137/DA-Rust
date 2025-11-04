use clap::{arg, Args, Parser, Subcommand};
use serde::{Serialize, Deserialize};

use crate::math::{EigenConfig, SplConfig, Jobz, Uplo};
use crate::{DEFAULT_STEP_NUM,
            DEFAULT_SYSTEM_SIZE,
            DEFAULT_INTERACTION_STRENGTH};

#[derive(Parser, Debug)]
#[command(
author="Thomas Prosser <thomas@prosser.ch>",
version,
about = "Simulation of kicked Bose-Einstein Condensates",
long_about = None)]
#[command(subcommand_precedence_over_arg = true)]
pub struct Cli {

    #[arg(short, long, global=true)]
    pub config: Option<String>,

    #[command(flatten)]
    pub global: GlobalConfig,

    #[command(subcommand)]
    pub alg: AlgorithmSubcommand,

}

/// Learned from ChatGPT: Seems Serde does not look at Clap
/// So default values decoupled from serializing/deserializing
///
/// What AI did not mention though is that Serde needs Strings
/// or functions for its default argument
fn default_step_num() -> usize { DEFAULT_STEP_NUM }
fn default_system_size() -> f64 { DEFAULT_SYSTEM_SIZE }
fn default_interaction_strength() -> f64 { DEFAULT_INTERACTION_STRENGTH }

#[derive(Args, Debug, Clone, Serialize, Deserialize, Default)]
pub struct GlobalConfig {

    /// Optional output file name (required if flag is used)
    #[arg(short, long, requires="format", value_name = "FILE")]
    pub output: Option<String>,

    /// File format. Only mandatory if output is set
    #[arg(short, long, value_parser=["png", "svg"])]
    pub format: Option<String>,

    /// Number of steps to use
    #[arg(short='n', long, default_value_t=DEFAULT_STEP_NUM, global=true)]
    #[serde(default = "default_step_num")]
    pub step_num: usize,

    /// Length of the system
    #[arg(short, long, default_value_t=DEFAULT_SYSTEM_SIZE, global=true)]
    #[serde(default = "default_system_size")]
    pub system_size: f64,

    /// Interaction strength of the bosons in the BEC, defaults to 0.0
    #[arg(short = 'g', long, default_value_t=DEFAULT_INTERACTION_STRENGTH, global=true)]
    #[serde(default = "default_interaction_strength")]
    pub interaction_strength: f64,

    /// Puts the BEC into an optical trap, defaults to false
    #[arg(short, long, default_value_t=false, global=true)]
    pub trap: bool,

    /// Adds a harmonic optical lattice, defaults to false
    #[arg(short, long, default_value_t=false, global=true)]
    pub lattice: bool,

}

#[derive(Subcommand, Debug, Serialize, Deserialize)]
pub enum AlgorithmSubcommand {
    Eig(EigenConfigArgs),
    Spl(SplConfigArgs),
}

#[derive(Args, Debug, Serialize, Deserialize)]
pub struct EigenConfigArgs {

    #[arg(long, value_enum, default_value_t = Jobz::WithEigenvectors)]
    pub mode: Jobz,

    #[arg(long, value_enum, default_value_t = Uplo::UpperTriangle)]
    pub symmetry: Uplo,

}

#[derive(Args, Debug, Serialize, Deserialize)]
pub struct SplConfigArgs {

    #[arg(long, default_value_t = 0.1)]
    pub dt: f64,

    #[arg(short, long, default_value_t = 1.0)]
    pub omega: f64,

    #[arg(short, long)]
    pub imag_time: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum AlgorithmConfig {
    Eig(EigenConfig),
    Spl(SplConfig),
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FullConfig {
    pub global: GlobalConfig,
    pub algorithm: AlgorithmConfig,
}
