 use num::pow;
 use crate::GlobalConfig;

 pub const PI: f64 = std::f64::consts::PI;
 pub const FRAC_ROOT_TWO_PI: f64 = 0.398942280401432677939946059934381868_f64;

#[derive(Clone, Debug)]
pub struct Hamiltonian {
    pub operator: Vec<f64>,
    pub interaction_strength: f64,
    pub trap: bool,
    pub lattice: bool,
}
impl Hamiltonian {
    pub fn new(config: &GlobalConfig, interaction_strength: f64, trap: bool, lattice: bool) -> Hamiltonian {
        let config = config.clone();
        let num_steps = config.step_num as usize;
        let fnum_steps= num_steps as f64;
        let system_width = config.clone().system_size;
        let operator = init_operator(&num_steps, &system_width, &fnum_steps, &interaction_strength, &trap, &lattice, |row, col| {
            let step_size = system_width / fnum_steps;
            if row == col {
                let wave_number = ( PI * 40. ) / system_width.clone();
                let xpos = position(&col, &system_width, &fnum_steps);
                let val = 1./(&step_size * &step_size)
                    + interaction_strength * pow(FRAC_ROOT_TWO_PI * f64::exp(- pow(xpos, 2) / 2.), 2)
                    + potential(&xpos, &wave_number, lattice, trap);
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
    pub diag: Option<Vec<f64>>,
    pub offdiag: Option<Vec<f64>>,
}

#[derive(Clone, Debug)]
pub struct TridiagHamiltonian {
    pub vectors: TridiagOperator,
    pub interaction_strength: f64,
    pub trap: bool,
    pub lattice: bool,
}
impl TridiagHamiltonian {
    pub fn new(config: &GlobalConfig) -> TridiagHamiltonian {
        let config = config.clone();
        let n = config.step_num;
        let system_width = config.system_size;
        let interaction_strength = config.interaction_strength;
        let trap = config.trap;
        let lattice = config.lattice;
        let fnum_steps= n as f64;
        let step_size = system_width / fnum_steps;

        let d: Vec<f64> = vec![0.0; n]
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let wave_number = ( PI * 40. ) / system_width.clone();
                let xpos = position(&idx, &system_width, &fnum_steps);
                1./(&step_size * &step_size)
                    + interaction_strength.clone() * pow(FRAC_ROOT_TWO_PI * f64::exp(- pow(xpos, 2) / 2.), 2)
                    + potential(&xpos, &wave_number, trap, lattice)

            })
            .collect();

        let e:Vec<f64> = vec![0.0; &n-1]
            .iter()
            .map(|_|{
                -0.5 * 1./(pow(step_size, 2))
            }).collect();

        let diag = if d.is_empty() { None } else { Some(d)};
        let offdiag = if e.is_empty() { None } else { Some(e)};

        let vectors = TridiagOperator{diag, offdiag };

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

pub fn potential(location: &f64, wave_number: &f64, trap: bool, lattice: bool) -> f64 {
    let sinx = f64::sin(wave_number * location);

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

 #[cfg(test)]
 mod tests {
     use crate::physics::hamiltonians::{position, potential};

     #[test]
     fn test_position() {
         let pos = position(&5, &10., &10.);
         let pos2 = position(&10, &10., &10.);

         assert_eq!(pos, 0.0);
         assert_eq!(pos2, -5.0);
     }

     // Trap only, Harmonic Oscillator
     #[test]
     fn test_potential() {
         let osc_pot = potential(&0., &1., &true, &false);
         let osc_pot2 = potential(&-5.0, &1., &true, &false);

         assert_eq!(osc_pot, 0.);
         assert_eq!(osc_pot2, 12.5);
     }

     // Lattice only. Sin^2 -like potential
     #[test]
     fn test_lat_potential() {
         let lat_pot = potential(&0., &1., &false, &true);

         assert_eq!(lat_pot, 0.);
     }
 }
