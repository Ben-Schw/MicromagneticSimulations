#ifndef FIELD_H
#define FIELD_H

#include "vector.h"
#include "parameters.h"
#include <random>

// Anisotropy field for uniaxial anisotropy
static inline Vec3 anisotropy_field_uniax(const Vec3& m) {
    Vec3 u = llg::anis_u.normalized();
    double mdotu = m.dot(u);

    // B_ani = (2 Ku / mu) (m·u) u
    const double pref = (2.0 * llg::Ku) / llg::mu;
    return u * (pref * mdotu);
}

// Thermal field for sLLG (Gaussian, independent components)
static inline Vec3 thermal_field(std::mt19937& rng)
{
    // sigma = sqrt( 2 α kB T / (γ μ dt) )
    const double num = 2.0 * llg::alpha * llg::kB * llg::T;
    const double den = llg::gamma * llg::mu * llg::dt;
    const double sigma = std::sqrt(num / den);

    std::normal_distribution<double> N01(0.0, 1.0);
    return Vec3(sigma * N01(rng), sigma * N01(rng), sigma * N01(rng));
}

// Get external magnetic field at time t
static inline Vec3 get_magnetic_field(double t) {
    constexpr double pi = 3.14159265358979323846;
    const double omega = 2.0 * pi * llg::f_B;
    const double arg = omega * t + llg::phase;

    switch (llg::field_mode) {
        case llg::FieldType::Constant:
            return llg::B0;

        case llg::FieldType::Oscillating:
            return llg::B0 + llg::B_var * std::sin(arg);

        // different planes of rotation
        case llg::FieldType::Rotating_XY:
            return llg::B0 + Vec3(
                llg::B_var.x * std::cos(arg),
                llg::B_var.y * std::sin(arg),
                0.0
            );
        case llg::FieldType::Rotating_YZ:
            return llg::B0 + Vec3(
                0.0,
                llg::B_var.y * std::cos(arg),
                llg::B_var.z * std::sin(arg)
            );
        case llg::FieldType::Rotating_XZ:
            return llg::B0 + Vec3(
                llg::B_var.x * std::cos(arg),
                0.0,
                llg::B_var.z * std::sin(arg)
            );
        case llg::FieldType::Pulse:
            if (t >= llg::time_delay && t <= llg::time_delay + llg::pulse_duration) {
                return llg::B0 + llg::B_var;
            } else {
                return llg::B0;
            }

        default:
            return llg::B0; // fallback
    }
}

// Get total effective field at time t for magnetization m
static inline Vec3 get_effective_field(double t, const Vec3& m) {
    const Vec3 Bext = get_magnetic_field(t);

    if (llg::use_anisotropy) {
        Bext += anisotropy_field_uniax(m);
    }

    return Bext;
}

#endif