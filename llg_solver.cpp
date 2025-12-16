#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

#include "vector.h"
#include "parameters.h"
#include "magnetic_field.h"

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

static inline double thermal_sigma(double alpha, double gamma,
                                      double T, double mu, double dt)
{
    // sigma = sqrt( 2 α kB T / (γ μ dt) )
    const double num = 2.0 * alpha * llg::kB * T;
    const double den = gamma * mu * dt;
    return std::sqrt(num / den);
}

// One Heun step for the stochastic LLG
static inline Vec3 heun_step_sllg(const Vec3& m, double t, double dt, const Vec3& Bth,
                                  double gamma, double alpha)
{
    // field at time t
    Vec3 B0 = get_magnetic_field(t);

    // predictor
    Vec3 f0 = llg_time_derivative(m, B0 + Bth, gamma, alpha);
    Vec3 m_pred = (m + f0 * dt).normalized();

    // magnetic field at t + dt
    Vec3 B1 = get_magnetic_field(t + dt);

    // corrector
    Vec3 f1 = llg_time_derivative(m_pred, B1 + Bth, gamma, alpha);

    Vec3 m_next = m + (f0 + f1) * (0.5 * dt);
    return m_next.normalized();
}

int main() {
    std::mt19937 rng(llg::rng_seed);
    std::normal_distribution<double> N01(0.0, 1.0);

    const double sigma = thermal_sigma(llg::alpha, llg::gamma,
                                       llg::T, llg::mu, llg::dt);

    Vec3 m = llg::m_init.normalized();

    std::ofstream out("output_files/traj.csv");
    out << "t,mx,my,mz,Bx,By,Bz\n";
    out << std::setprecision(12);

    for (int i = 0; i < llg::nsteps; ++i) {
        double t = i * llg::dt;

        Vec3 B = get_magnetic_field(t);

        // Gaussian thermal field (Tesla)
        Vec3 Bth(sigma * N01(rng), sigma * N01(rng), sigma * N01(rng));

        out << t << "," << m.x << "," << m.y << "," << m.z << ","
            << B.x << "," << B.y << "," << B.z << "\n";

        m = heun_step_sllg(m, t, llg::dt, Bth, llg::gamma, llg::alpha);
    }
    return 0;
}
