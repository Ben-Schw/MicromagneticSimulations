#ifndef LLG_PARAMS_H
#define LLG_PARAMS_H

#include <cmath>
#include "vector.h"

namespace llg {
    // simulation parameters
    inline constexpr double gamma = 1.760859e11; // rad/(s*T)
    inline constexpr double alpha = 0.05;

    inline constexpr double dt = 1e-16; // s
    inline constexpr int    nsteps = 800000;

    // ---- Lattice ----
    inline constexpr int Nx = 16;
    inline constexpr int Ny = 16;
    inline constexpr int Nz = 16;
    inline constexpr bool periodic = true;

    // Exchange coupling (effective-field scale):
    // J > 0 : Ferromagnet
    // J < 0 : Antiferromagnet
    inline constexpr double J = +0.2;

    // Init
    enum class InitType { Uniform, AFM_Checkerboard, Random };
    inline constexpr InitType init_type = InitType::AFM_Checkerboard;

    inline const Vec3 m_init{1.0, 0.0, 0.0}; // used for Uniform init

    // ---- Output ----
    enum class OutputAxis { X, Y, Z };
    inline constexpr OutputAxis output_axis = OutputAxis::Z;

    // external field (Tesla)
    enum class FieldType { Constant, Oscillating, Rotating_XY, Rotating_YZ, Rotating_XZ, Pulse }; // rotation in respective plane anticlockwise when looking along positive axis

    inline const Vec3 B_var{0.0, 0.0, 0.0};
    inline constexpr double f_B = 1e13; // Hz
    inline constexpr double phase = 0.0; // rad, phase = 0: sinusoidal start, phase = pi/2: cosinusoidal start
    inline const Vec3 B0{0.0, 0.0, 15.0};
    inline const double pulse_duration = 5e-12; // only for FieldType::Pulse
    inline const double time_delay = 10e-12;      // only for FieldType::Pulse

    inline constexpr FieldType field_mode = FieldType::Constant;

    // thermal
    inline constexpr double kB  = 1.380649e-23;      // J/K
    inline constexpr double muB = 9.2740100783e-24;  // J/T

    inline constexpr double T = 2.0; // K

    // spin parameters
    inline constexpr double g = 2.0;
    inline constexpr double S = 1.0;               // spin quantum number
    inline constexpr double mu = g * muB * S;      // J/T

    inline constexpr unsigned int rng_seed = 12345;
}


#endif