use clap::{arg, Parser};
use crate::math::Jobz;

#[derive(Parser, Debug)]
#[command(
author="Thomas Prosser <thomas@prosser.ch>",
version,
about = "Simulation of kicked Bose-Einstein Condensates",
long_about = None)]
pub struct Cli {
    /// Optional output file name (required if flag is used)
    #[arg(short, long, requires="format", value_name = "FILE")]
    pub output: Option<String>,

    /// File format. Only mandatory if output is set
    #[arg(short, long, value_parser=["png", "svg"])]
    pub format: Option<String>,

    /// Number of steps to use
    #[arg(short='n', long)]
    pub step_num: usize,

    /// Length of the system
    #[arg(short, long)]
    pub system_size: f64,

    /// Interaction strength of the bosons in the BEC, defaults to 0.0
    #[arg(short = 'g', long, default_value_t=0.0)]
    pub interaction_strength: f64,

    /// Puts the BEC into an optical trap, defaults to false
    #[arg(short, long, default_value_t=false)]
    pub trap: bool,

    /// Adds a harmonic optical lattice, defaults to false
    #[arg(short, long, default_value_t=false)]
    pub lattice: bool,

    /// Defaults to include Eigenvectors
    #[arg(value_enum, value_name = "MODE", default_value_t = Jobz::WithEigenvectors)]
    pub mode: Jobz
}

