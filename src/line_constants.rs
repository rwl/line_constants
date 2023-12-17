// Copyright (C) 2017 Julius Susanto
// Copyright (C) 2023 Richard Lincoln
// All rights reserved.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::f64::consts::PI;

use num_complex::Complex;

const MU0: f64 = 4.0 * PI * 1e-7; // Permeability of free space

trait Sqr {
    fn sqr(&self) -> Self;
}

impl Sqr for f64 {
    fn sqr(&self) -> f64 {
        self.powi(2)
    }
}

pub enum ConductorType {
    Solid,
    Tubular,
}

/// Calculates internal inductance of solid or tubular conductor
/// Note that calculations assume uniform current distribution in the conductor, thus conductor stranding is not taken into account.
///
/// Usage:
///     L_int = calc_L_int(type, r, q)
///
/// where   type is 'solid' or 'tube'
///         r is the radius of the conductor (mm)
///         q is the radius of the inner tube (mm)
///
/// Returns:
///         L_int the internal inductance of the conductor (mH/m)
pub fn calc_l_int(conductor_type: ConductorType, r: f64, q: f64) -> f64 {
    match conductor_type {
        ConductorType::Solid => MU0 / 8.0 / PI,
        ConductorType::Tubular => {
            MU0 / 2.0 / PI
                * (q.powi(4) / (r.sqr() - q.sqr()).sqr() * (r / q).ln()
                    - (3.0 * q.sqr() - r.sqr()) / (4.0 * (r.sqr() - q.sqr())))
        }
    }
}

/// Calculates geometric mean radius (GMR) of solid or tubular conductor
/// Note that calculations assume uniform current distribution in the conductor,
/// thus conductor stranding is not taken into account.
///
/// Usage:
///     GMR = calc_GMR(type, r, q)
///
/// where   type is 'solid' or 'tube'
///         r is the radius of the conductor (mm)
///         q is the radius of the inner tube (mm)
///
/// Returns:
///         GMR the geometric mean radius (mm)
pub fn calc_gmr(conductor_type: ConductorType, r: f64, q: f64) -> f64 {
    match conductor_type {
        ConductorType::Solid => r * f64::exp(-0.25),
        ConductorType::Tubular => {
            r * f64::exp(
                (3.0 * q.sqr() - r.sqr()) / (4.0 * (r.sqr() - q.sqr()))
                    - q.powi(4) / (r.sqr() - q.sqr()).sqr() * f64::ln(r / q),
            )
        }
    }
}

pub enum ImpedanceType {
    Self_,
    Mutual,
}

/// Calculates Carson's earth return correction factors Rp and Xp for both self and mutual terms.
/// The number of terms evaluated in the infinite loop is based on convergence to the desired error tolerance.
///
/// Usage:
///     Rp, Xp = carsons(type, h_i, h_k, x_ik, f, rho, err_tol)
///
/// where   type is 'self' or 'mutual'
///         h_i is the height of conductor i above ground (m)
///         h_k is the height of conductor k above ground (m)
///         x_ik is the horizontal distance between conductors i and k (m)
///         f is the frequency (Hz)
///         rho is the earth resistivity (Ohm.m)
///         err_tol is the error tolerance for the calculation (default = 1e-6)
///
/// Returns:
///         Rp, Xp the Carson earth return correction factors (in Ohm/km)
pub fn carsons(
    term_type: ImpedanceType,
    h_i: f64,
    h_k: f64,
    x_ik: f64,
    f: f64,
    rho: f64,
    err_tol: Option<f64>, /* =1e-6 */
) -> (f64, f64) {
    // Geometrical calculations
    let (d, cos_phi, sin_phi, phi) = match term_type {
        ImpedanceType::Self_ => {
            //     if type == 'self':
            let d = 2.0 * h_i;
            let cos_phi = 1.0;
            let sin_phi = 0.0;
            let phi = 0.0;
            (d, cos_phi, sin_phi, phi)
        }
        ImpedanceType::Mutual => {
            let d = f64::sqrt((h_i + h_k).sqr() + x_ik.sqr()); // Distance between conductor i and image of conductor k
            let cos_phi = (h_i + h_k) / d;
            let sin_phi = (x_ik) / d;
            // let phi = np.arccos(cos_phi);
            let phi = f64::acos(cos_phi);
            (d, cos_phi, sin_phi, phi)
        }
    };

    // Initialise parameters
    let mut i = 1.0;
    let mut err = 1.0;
    let mut sgn = 1.0;

    // Initial values and constants for calculation
    let omega = 2.0 * PI * f;
    let a = 4.0 * PI * f64::sqrt(5.0) * 1e-4 * d * f64::sqrt(f / rho);
    let mut acosphi = a * cos_phi;
    let mut asinphi = a * sin_phi;
    let mut b = vec![f64::sqrt(2.0) / 6.0, 1.0 / 16.0];
    let mut c = vec![0.0, 1.3659315];
    let mut d = b.iter().map(|&bi| PI / 4.0 * bi).collect::<Vec<f64>>();

    // First two terms of carson correction factor
    let mut rp = PI / 8.0 - b[0] * acosphi;
    let mut xp = 0.5 * (0.6159315 - f64::ln(a)) + b[0] * acosphi;

    // Loop through carson coefficient terms starting with i = 2
    while err > err_tol.unwrap_or(1e-6) {
        // let term = f64::mod(i, 4);
        let term = i % 4.0;
        // Check sign for b term
        if term == 0.0 {
            sgn = -1.0 * sgn;
        }
        // Calculate coefficients
        let bi = b[i as usize - 1] * sgn / ((i + 1.0) * (i + 3.0));
        let ci = c[i as usize - 1] + 1.0 / (i + 1.0) + 1.0 / (i + 3.0);
        let di = PI / 4.0 * bi;
        b.push(bi);
        c.push(ci);
        d.push(di);

        // Recursively calculate powers of acosphi and asinphi
        let acosphi_prev = acosphi;
        let asinphi_prev = asinphi;
        acosphi = (acosphi_prev * cos_phi - asinphi_prev * sin_phi) * a;
        asinphi = (acosphi_prev * sin_phi + asinphi_prev * cos_phi) * a;

        let rp_prev = rp;
        let xp_prev = xp;

        // First term
        if term == 0.0 {
            rp = rp - bi * acosphi;
            xp = xp + bi * acosphi;
        }
        // Second term
        else if term == 1.0 {
            rp = rp + bi * ((ci - a.ln()) * acosphi + phi * asinphi);
            xp = xp - di * acosphi;
        }
        // Third term
        else if term == 1.0 {
            rp = rp + bi * acosphi;
            xp = xp + bi * acosphi;
        }
        // Fourth term
        else {
            rp = rp - di * acosphi;
            xp = xp - bi * ((ci - a.ln()) * acosphi + phi * asinphi);
        }

        // i = i = 1;
        i = i + 1.0;
        err = f64::sqrt((rp - rp_prev).sqr() + (xp - xp_prev).sqr());
    }

    rp = 4.0 * omega * 1e-04 * rp;
    xp = 4.0 * omega * 1e-04 * xp;

    return (rp, xp);
}

/// Calculates self impedance term (in Ohm/km)
/// NOTE: No allowance has been made for skin effects
///
/// Usage:
///     self_Z = calc_self_Z(R_int, cond_type, r, q, h_i, f, rho, err_tol=1e-6)
///
/// where   R_int is the AC conductor resistance (Ohm/km)
///         cond_type is the conductor type ('solid' or 'tube')
///         r is the radius of the conductor (mm)
///         q is the radius of the inner tube (mm)
///         h_i is the height of conductor i above ground (m)
///         f is the frequency (Hz)
///         rho is the earth resistivity (Ohm.m)
///         err_tol is the error tolerance for the calculation (default = 1e-6)
///
/// Returns:
///         self_Z the self impedance term of line impedance matrix (Ohm/km)
pub fn calc_self_z(
    r_int: f64,
    cond_type: ConductorType,
    r: f64,
    q: f64,
    h_i: f64,
    f: f64,
    rho: f64,
    err_tol: Option<f64>,
) -> Complex<f64> {
    // Constants
    let omega = 2.0 * PI * f; // Nominal angular frequency

    // Calculate internal conductor reactance (in Ohm/km)
    let x_int = 1000.0 * omega * calc_l_int(cond_type, r, q);

    // Calculate geometrical reactance (in Ohm/km)
    let x_geo = 1000.0 * omega * MU0 / 2.0 / PI * (2.0 * h_i / r * 1000.0).ln();

    // Calculate Carson's correction factors (in Ohm/km)
    let (rp, xp) = carsons(ImpedanceType::Self_, h_i, 0.0, 0.0, f, rho, err_tol);

    Complex::new(r_int + rp, x_int + x_geo + xp)
}

/// Calculates mutual impedance term (in Ohm/km)
///
/// Usage:
///     mutual_Z = calc_mutual_Z(cond_type, r, q, h_i, h_k, x_ik, f, rho, err_tol=1e-6)
///
/// where   cond_type is the conductor type ('solid' or 'tube')
///         r is the radius of the conductor (mm)
///         q is the radius of the inner tube (mm)
///         h_i is the height of conductor i above ground (m)
///         h_k is the height of conductor k above ground (m)
///         x_ik is the horizontal distance between conductors i and k (m)
///         f is the frequency (Hz)
///         rho is the earth resistivity (Ohm.m)
///         err_tol is the error tolerance for the calculation (default = 1e-6)
///
/// Returns:
///         mutual_Z the self impedance term of line impedance matrix (Ohm/km)
pub fn calc_mutual_z(
    _cond_type: ConductorType,
    _r: f64,
    _q: f64,
    h_i: f64,
    h_k: f64,
    x_ik: f64,
    f: f64,
    rho: f64,
    err_tol: Option<f64>,
) -> Complex<f64> {
    // Constants
    let omega = 2.0 * PI * f; // Nominal angular frequency
    let d_img = f64::sqrt((h_i + h_k).sqr() + x_ik.sqr()); // Distance between conductor i and image of conductor k
    let d = f64::sqrt((h_i - h_k).sqr() + x_ik.sqr()); // Distance between conductors i and k

    // Calculate geometrical mutual reactance (in Ohm/km)
    let x_geo = 1000.0 * omega * MU0 / 2.0 / PI * f64::ln(d_img / d);

    // Calculate Carson's correction factors (in Ohm/km)
    let (r_p, x_p) = carsons(ImpedanceType::Mutual, h_i, h_k, x_ik, f, rho, err_tol);

    Complex::new(r_p, x_geo + x_p)
}
