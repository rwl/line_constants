//! Overhead Line Constants Calculation Library
//!
//! Author: Julius Susanto  
//! Last update: January 2017
//!
//! Reference:
//!  1. Dommel, H. W., "Electromagnetic Transients Program Reference Manual (EMTP Theory Book)", Chapter 4, "Overhead Transmission Lines".
//!  2. Arrillaga, J., and Watson, N. R., "Computer Modelling of Electrical Power Systems", 2nd Edition, Wiley, 2005, Chapter 2.6
//!
//! Functions:
//!  - `calc_l_int`          Calculates internal inductance of solid or tubular conductor
//!  - `calc_gmr`            Calculates geometric mean radius (GMR) of solid or tubular conductor
//!  - `carsons`             Calculates Carson's earth return correction factors Rp and Xp for self or mutual terms
//!  - `calc_self_z`         Calculates self impedance term (in Ohm/km)
//!  - `calc_mutual_z`       Calculates mutual impedance term (in Ohm/km)
//!  - `calc_dubanton_z`     Calculates Dubanton approximation for self or mutual impedance (in Ohm/km)
//!  - `calc_z_matrix`       Calculates primitive impedance matrix
//!  - `calc_y_matrix`       Calculates primitive admittance matrix
//!  - `calc_kron_z`         Calculates Kron reduced matrix
//!
//! Features currently not supported:
//!  - Skin effect calculation for internal self impedance
//!  - Input data validation, e.g. check that all phase conductor vectors are the same size
mod line_constants;

pub use line_constants::*;
