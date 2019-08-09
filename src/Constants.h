#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <random>

// Global Constants
const double k_b = 1.38e-23;
const double T = 87;
const double M_a = 6.6e-26;
const double m_e = 9.1e-31;
const double e = 1.6e-19;

// Parameters
const bool is_gas = true;
const bool use_recursive_ionization = false;
const bool uniform_field = true;
const bool interactions = true;

const bool single_tip = true;
const double spawn_height_scale = 0.6;
const double spawn_width_scale = 0.1;

// Liquid number density
// const double n = 2.11e22;

// Gas number densities
// const double n = 1.15e19; // at 2 PSI (103 Torr)
// const double n = 4.06e19; // at 7 PSI (362 Torr) (really 4.06 cm^-3) (was 4.02)
const double n = 8.61e19; // at 15 PSI (776 Torr) (really 8.81 cm^-3)

#endif