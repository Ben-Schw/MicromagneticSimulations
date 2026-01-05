/*
  LLG_solver_spin_grid.cpp

  Solve the stochastic Landau-Lifshitz-Gilbert equation for a grid of spins.

  Uses Heun's method (predictor-corrector) for integration.
*/
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

#include "vector.h"
#include "parameters.h"
#include "effective_field.h"
#include "grid.h"


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
static inline Vec3 heun_step_sllg(SpinGrid& grid, int x, int y, int z, Vec3& m, double t, double dt,
                                  double gamma, double alpha, std::mt19937& rng)
{
    // thermal field (same for predictor and corrector)
    Vec3 Bth{0.0, 0.0, 0.0};
    if (llg::use_thermal_field) {
        Vec3 Bth = thermal_field(rng);
    }

    // field at time t (deterministic)
    const Vec3 B0 = get_effective_field(t, m);

    // exchange field contribution for multi-spin systems
    const Vec3 B_exch_0 = exchange_field(grid, x, y, z, llg::mu);


    // predictor
    const Vec3 f0 = llg_time_derivative(m, B0 + Bth + B_exch_0, gamma, alpha);
    const Vec3 m_pred = (m + f0 * dt).normalized();

    // field at time t + dt (deterministic, depends on predicted state)
    const Vec3 B1 = get_effective_field(t + dt, m_pred);
    const Vec3 B_exch_1 = exchange_field(grid, x, y, z, llg::mu);


    // corrector (IMPORTANT: use B1 here!)
    const Vec3 f1 = llg_time_derivative(m_pred, B1 + Bth + B_exch_1, gamma, alpha);

    const Vec3 m_next = m + (f0 + f1) * (0.5 * dt);
    return m_next.normalized();
}


int main() {
    std::mt19937 rng(llg::rng_seed);
    SpinGrid grid(llg::Nx, llg::Ny, llg::Nz, llg::bc);
    grid.initialize_all(llg::m_init);
    std::ofstream out("traj_spin_grid.csv");
    out << "t,x,y,z,mx,my,mz,Bx,By,Bz\n";
    out << std::setprecision(12);   
    for (int step = 0; step < llg::nsteps; ++step) {
        double t = step * llg::dt;

        for (int z = 0; z < llg::Nz; ++z) {
            for (int y = 0; y < llg::Ny; ++y) {
                for (int x = 0; x < llg::Nx; ++x) {
                    int idx = grid.index_raw(x, y, z);
                    Vec3 m = grid.spin_linear(idx);

                    Vec3 B = get_magnetic_field(t);

                    out << t << "," << x << "," << y << "," << z << ","
                        << m.x << "," << m.y << "," << m.z << ","
                        << B.x << "," << B.y << "," << B.z << "\n";

                    Vec3 m_next = heun_step_sllg(grid, x, y, z, m, t, llg::dt, llg::gamma, llg::alpha, rng);
                    grid.at_raw(x, y, z) = m_next;
                }
            }
        }
    }
    return 0;
}

