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

// -----------------------------------------------------------------------------
// Perform one Heun (predictor–corrector) step for the full spin grid.
//
// The update is fully synchronous:
//   - all effective fields in the predictor are computed from the old state
//   - all effective fields in the corrector are computed from the predictor
//
// Thermal noise (if enabled) is drawn once per spin and reused in both
// predictor and corrector stages, corresponding to a Stratonovich treatment.
// -----------------------------------------------------------------------------
static void heun_step_grid(SpinGrid& grid, double t, double dt, std::mt19937& rng)
{
    // Store the state at time t
    SpinGrid old = grid;

    // Predictor grid (same geometry and boundary conditions)
    SpinGrid pred(old.Nx, old.Ny, old.Nz, old.bc);

    const int N = old.Nx * old.Ny * old.Nz;

    // Thermal field storage (one realization per spin and timestep)
    std::vector<Vec3> Bth_store(static_cast<size_t>(N), Vec3(0.0, 0.0, 0.0));

    if (llg::use_thermal_field) {
        for (int idx = 0; idx < N; ++idx) {
            Bth_store[static_cast<size_t>(idx)] = thermal_field(rng);
        }
    }

    // ------------------------
    // Predictor (Euler step)
    // ------------------------
    for (int z = 0; z < old.Nz; ++z) {
        for (int y = 0; y < old.Ny; ++y) {
            for (int x = 0; x < old.Nx; ++x) {
                const int idx = old.index_raw(x, y, z);

                const Vec3 m   = old.spin_linear(idx);
                const Vec3 Bth = Bth_store[static_cast<size_t>(idx)];

                // Effective field contributions evaluated at time t
                const Vec3 B0   = get_effective_field(t, m);
                const Vec3 Bex0 = exchange_field(old, x, y, z, llg::mu);

                // Time derivative dm/dt from the LLG equation
                const Vec3 f0 = llg_time_derivative(
                    m,
                    B0 + Bex0 + Bth,
                    llg::gamma,
                    llg::alpha
                );

                // Predictor state (normalized explicitly)
                pred.at_raw(x, y, z) = (m + f0 * dt).normalized();
            }
        }
    }

    // ------------------------
    // Corrector (Heun step)
    // ------------------------
    for (int z = 0; z < old.Nz; ++z) {
        for (int y = 0; y < old.Ny; ++y) {
            for (int x = 0; x < old.Nx; ++x) {
                const int idx = old.index_raw(x, y, z);

                const Vec3 m      = old.spin_linear(idx);
                const Vec3 m_pred = pred.spin_linear(idx);
                const Vec3 Bth    = Bth_store[static_cast<size_t>(idx)];

                // Effective fields at t and t + dt
                const Vec3 B0   = get_effective_field(t, m);
                const Vec3 B1   = get_effective_field(t + dt, m_pred);

                const Vec3 Bex0 = exchange_field(old,  x, y, z, llg::mu);
                const Vec3 Bex1 = exchange_field(pred, x, y, z, llg::mu);

                const Vec3 f0 = llg_time_derivative(
                    m,
                    B0 + Bex0 + Bth,
                    llg::gamma,
                    llg::alpha
                );

                const Vec3 f1 = llg_time_derivative(
                    m_pred,
                    B1 + Bex1 + Bth,
                    llg::gamma,
                    llg::alpha
                );

                // Final update (normalized explicitly)
                grid.at_raw(x, y, z) = (m + (f0 + f1) * (0.5 * dt)).normalized();
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Compute the average magnetization of a layer perpendicular to the given axis.
// The layer index refers to the coordinate along the chosen axis.
// -----------------------------------------------------------------------------
static Vec3 layer_magnetization(const SpinGrid& grid, llg::OutputAxis axis, int layer)
{
    Vec3 sum(0.0, 0.0, 0.0);
    int count = 0;

    switch (axis) {
        case llg::OutputAxis::X:
            for (int z = 0; z < grid.Nz; ++z) {
                for (int y = 0; y < grid.Ny; ++y) {
                    sum += grid.at_raw(layer, y, z);
                    ++count;
                }
            }
            break;

        case llg::OutputAxis::Y:
            for (int z = 0; z < grid.Nz; ++z) {
                for (int x = 0; x < grid.Nx; ++x) {
                    sum += grid.at_raw(x, layer, z);
                    ++count;
                }
            }
            break;

        case llg::OutputAxis::Z:
            for (int y = 0; y < grid.Ny; ++y) {
                for (int x = 0; x < grid.Nx; ++x) {
                    sum += grid.at_raw(x, y, layer);
                    ++count;
                }
            }
            break;
    }

    return sum * (1.0 / static_cast<double>(count));
}

int main()
{
    // Random number generator for thermal noise
    std::random_device rd;
    std::mt19937 rng(rd());

    // Initialize spin grid
    SpinGrid grid(llg::Nx, llg::Ny, llg::Nz, llg::bc);

    grid.initialize_all(llg::m_init);

    std::ofstream out("spin_grid_output.csv");
    out << "t,T,z,mx,my,mz,Bx,By,Bz\n";

    // Time integration loop
    for (int step = 0; step < llg::nsteps; ++step) {
        const double t = step * llg::dt;
        const Vec3 B = get_magnetic_field(t);
        const double temp = llg::T;

        // Output layer-resolved magnetization if requested
        const int n_layers =
            (llg::output_axis == llg::OutputAxis::X) ? grid.Nx :
            (llg::output_axis == llg::OutputAxis::Y) ? grid.Ny :
                                    grid.Nz;
                                    
        for (int i = 0; i < n_layers; ++i) {
            Vec3 m = layer_magnetization(grid, llg::output_axis, i);
            out << t << "," << temp << "," << i << ","
                << m.x << "," << m.y << "," << m.z << ","
                << B.x << "," << B.y << "," << B.z << "\n";
        }

        // Advance the full system by one time step
        heun_step_grid(grid, t, llg::dt, rng);

    }

    return 0;
}