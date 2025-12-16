#ifndef FIELD_H
#define FIELD_H

#include "vector.h"
#include "parameters.h"


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

#endif