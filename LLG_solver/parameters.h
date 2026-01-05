#ifndef LLG_PARAMS_H
#define LLG_PARAMS_H

#include <cmath>
#include <array>
#include "vector.h"

namespace llg {
    // simulation parameters
    inline constexpr double gamma = 1.760859e11; // rad/(s*T)
    inline constexpr double alpha = 0.05;    // dimensionless damping parameter
    inline constexpr double meV_to_J = 1.602176634e-22; // 1 meV in Joule


    inline constexpr double dt = 1e-16; // s
    inline constexpr int    nsteps = 100000;

    // --- Anisotropy ---
    inline constexpr bool   use_anisotropy = true;

    // Easy-axis direction (will be normalized in code, but keep it nice)
    inline const Vec3 anis_u{0.0, 0.0, 1.0};

    // Uniaxial anisotropy strength Ku in meV
    inline constexpr double Ku = 0.5 * meV_to_J; // meV

    // Exchange coupling (effective-field scale):
    // J_ij in meV between neighboring spins at (i,j,k) and (i+dx,j+dy,k+dz)
    struct Interaction {
        int dx, dy, dz;
        double J;
    };
    // Nearest-neighbor exchange interactions (J in meV)
    inline constexpr std::array<Interaction, 6> exchange_interactions {{
        {+1, 0, 0, 50.0 * meV_to_J},
        {-1, 0, 0, 50.0 * meV_to_J},
        {0, +1, 0, 50.0 * meV_to_J},
        {0, -1, 0, 50.0 * meV_to_J},
        {0, 0, +1, 50.0 * meV_to_J},
        {0, 0, -1, 50.0 * meV_to_J},
    }};

    // external field (Tesla)
    inline constexpr bool use_external_field = true;
    enum class FieldType { Constant, Oscillating, Rotating_XY, Rotating_YZ, Rotating_XZ, Pulse }; // rotation in respective plane anticlockwise when looking along positive axis

    inline const Vec3 B_var{0.0, 0.0, 0.0};
    inline constexpr double f_B = 1e13; // Hz
    inline constexpr double phase = 0.0; // rad, phase = 0: sinusoidal start, phase = pi/2: cosinusoidal start
    inline const Vec3 B0{0.0, 0.0, 0.0};
    inline const double pulse_duration = 5e-12; // only for FieldType::Pulse
    inline const double time_delay = 10e-12;      // only for FieldType::Pulse

    inline constexpr FieldType field_mode = FieldType::Constant;

    
    // thermal field
    inline constexpr bool use_thermal_field = true;
    inline constexpr double kB  = 1.380649e-23;      // J/K
    inline constexpr double muB = 9.2740100783e-24;  // J/T

    inline constexpr double T = 10.0; // K

    // spin parameters
    inline constexpr double g = 2.0;
    inline constexpr double S = 1.0;               // spin quantum number
    inline constexpr double mu = g * muB * S;      // J/T

    // ---- Lattice ----
    inline constexpr int Nx = 16;
    inline constexpr int Ny = 16;
    inline constexpr int Nz = 16;

    // Boundary conditions
    inline constexpr Vec3i bc{1, 1, 1}; // 0: open, 1: periodic

    // Init
    enum class InitType { Uniform, AFM_Checkerboard, Random };
    inline constexpr InitType init_type = InitType::AFM_Checkerboard;

    inline const Vec3 m_init{0.0, 0.0, 1.0}; // used for Uniform init

    // ---- Output ----
    enum class OutputAxis { X, Y, Z };
    inline constexpr OutputAxis output_axis = OutputAxis::Z;

    inline constexpr unsigned int rng_seed = 12345;
}


#endif