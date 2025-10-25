
 use num::pow;
 use crate::{FRAC_ROOT_TWO_PI, linalg};
 use linalg::EigenConfig;

#[derive(Clone, Debug)]
pub struct Hamiltonian {
    pub operator: Vec<f64>,
    pub interaction_strength: f64,
    pub trap: bool,
    pub lattice: bool,
}
impl Hamiltonian {
    pub fn new(config: &EigenConfig, interaction_strength: f64, trap: bool, lattice: bool) -> Hamiltonian {
        let config = config.clone();
        let num_steps = config.n as usize;
        let fnum_steps= num_steps as f64;
        let system_width = config.system_width;
        let operator = init_operator(&num_steps, &system_width, &fnum_steps, &interaction_strength, &trap, &lattice, |row, col| {
            let step_size = system_width / fnum_steps;
            if row == col {
                let wave_number = ( crate::PI * 40. ) / system_width.clone();
                let xpos = position(&col, &system_width, &fnum_steps);
                let val = 1./(&step_size * &step_size)
                    + interaction_strength * pow(FRAC_ROOT_TWO_PI * f64::exp(- pow(xpos, 2) / 2.), 2)
                    + potential(&xpos, &wave_number, &lattice, &trap);
                return val
            } else if (col as isize - row as isize).abs() == 1  {
                return -0.5 * 1./(pow(step_size, 2))
            }
            0.0
        });

        Hamiltonian{ operator, interaction_strength, trap, lattice }
    }



} // Impl Hamiltonian

#[derive(Clone, Debug)]
pub struct TridiagOperator {
    pub diag: Vec<f64>,
    pub offdiag: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct TridiagHamiltonian {
    pub vectors: TridiagOperator,
    pub interaction_strength: f64,
    pub trap: bool,
    pub lattice: bool,
}
impl TridiagHamiltonian {
    pub fn new(config: &EigenConfig, interaction_strength: f64, trap: bool, lattice: bool) -> TridiagHamiltonian {
        let config = config.clone();
        let n = config.n as usize;
        let system_width = config.system_width;
        let fnum_steps= n as f64;
        let step_size = system_width / fnum_steps;

        let d = vec![0.0; n]
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let wave_number = ( crate::PI * 40. ) / system_width.clone();
                let xpos = position(&idx, &system_width, &fnum_steps);
                1./(&step_size * &step_size)
                    + interaction_strength.clone() * pow(FRAC_ROOT_TWO_PI * f64::exp(- pow(xpos, 2) / 2.), 2)
                    + potential(&xpos, &wave_number, &lattice, &trap)

            })
            .collect();

        let e = vec![0.0; &n-1]
            .iter()
            .map(|_|{
                -0.5 * 1./(pow(step_size, 2))
            }).collect();

        let vectors = TridiagOperator{diag: d, offdiag: e };

        TridiagHamiltonian{ vectors, interaction_strength, trap, lattice}
    }
}


pub fn init_operator(num_steps: &usize, _system_width: &f64, _fnum_steps: &f64, _interaction_strength: &f64, _lattice: &bool, _trap: &bool, f: impl Fn(usize, usize) -> f64) -> Vec<f64> {
    let mut operator: Vec<f64> = Vec::new();
    for row in 0..*num_steps {
        for col in 0..*num_steps {
            let row = row.clone();
            operator.push( f(row, col));
        }
    };
    operator
}

pub fn position(index: &usize, system_width: &f64, fnum_steps: &f64) -> f64 {
    let idx = (*index).clone();
    (system_width * 0.5) - (idx as f64) * (system_width / fnum_steps)
}

pub fn potential(location: &f64, wave_number: &f64, trap: &bool, lattice: &bool) -> f64 {
    let sinx = f64::sin( wave_number * location );

    match (trap, lattice) {
        (true, false) => location * location * 0.5,
        (false, true) => {
            let pot = 0.5 * pow(sinx, 2) * location * location;
            pot
        },
        (true, true) => {
            let sinx_sq = pow(sinx, 2);
            let pot = 0.5 * location * location + 0.5 * wave_number * wave_number * &sinx_sq;
            pot
        },
        _ => 0.
    }
}

