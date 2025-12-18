#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

#include "vector.h"
#include "parameters.h"
#include "effective_field.h"

// dm/dt according to LLG
static inline Vec3 llg_time_derivative(
    const Vec3& m,
    const Vec3& B_eff,
    double gamma,
    double alpha)
{
    Vec3 mxB   = m.cross(B_eff);
    Vec3 mxmxB = m.cross(mxB);

    double prefactor = -gamma / (1.0 + alpha*alpha);
    return (mxB + mxmxB * alpha) * prefactor;
}



// One Heun step for the stochastic LLG (Stratonovich-style: same noise in predictor/corrector)
static inline Vec3 heun_step_sllg(const Vec3& m, double t, double dt,
                                  double gamma, double alpha, std::mt19937& rng)
{
    // thermal field (same for predictor and corrector)
    Vec3 Bth{0.0, 0.0, 0.0};
    if (llg::use_thermal_field) {
        Vec3 Bth = thermal_field(rng);
    }

    // field at time t (deterministic)
    const Vec3 B0 = get_effective_field(t, m);

    // predictor
    const Vec3 f0 = llg_time_derivative(m, B0 + Bth, gamma, alpha);
    const Vec3 m_pred = (m + f0 * dt).normalized();

    // field at time t + dt (deterministic, depends on predicted state)
    const Vec3 B1 = get_effective_field(t + dt, m_pred);

    // corrector (IMPORTANT: use B1 here!)
    const Vec3 f1 = llg_time_derivative(m_pred, B1 + Bth, gamma, alpha);

    const Vec3 m_next = m + (f0 + f1) * (0.5 * dt);
    return m_next.normalized();
}

int main() {
    std::mt19937 rng(llg::rng_seed);
    Vec3 m = llg::m_init.normalized();

    std::ofstream out("traj_single_spin.csv");
    out << "t,mx,my,mz,Bx,By,Bz\n";
    out << std::setprecision(12);

    for (int i = 0; i < llg::nsteps; ++i) {
        double t = i * llg::dt;

        Vec3 B = get_magnetic_field(t);

        out << t << "," << m.x << "," << m.y << "," << m.z << ","
            << B.x << "," << B.y << "," << B.z << "\n";
        
        // advance one time step
        m = heun_step_sllg(m, t, llg::dt, llg::gamma, llg::alpha, rng);
    }
    return 0;
}
