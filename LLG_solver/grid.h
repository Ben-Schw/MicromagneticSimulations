#ifndef GRID_H
#define GRID_H

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "vector.h"
#include "parameters.h"

class SpinGrid {
public:
    // By default, boundary conditions come from parameters.h
    SpinGrid(int Nx_, int Ny_, int Nz_, Vec3i bc_ = llg::bc) 
        : Nx(Nx_), Ny(Ny_), Nz(Nz_), bc(bc_), 
        spins(static_cast<size_t>(Nx_) * static_cast<size_t>(Ny_) * static_cast<size_t>(Nz))
    {
        if (Nx <= 0 || Ny <= 0 || Nz <= 0) {
            throw std::invalid_argument("Grid dimensions must be positive integers.");
        }
    }

    Vec3& at_raw(int x, int y, int z) { return spins[index_raw(x, y, z)]; }
    const Vec3& at_raw(int x, int y, int z) const { return spins[index_raw(x, y, z)]; }

    Vec3& spin_linear(int i) { return spins[static_cast<size_t>(i)]; }
    const Vec3& spin_linear(int i) const { return spins[static_cast<size_t>(i)]; }
    
    // Returns -1 if an open boundary is hit.
    int neighbor_index(int x, int y, int z, int dx, int dy, int dz) const
    {
        const int nx = apply_bc(x + dx, Nx, bc.periodic_x());
        const int ny = apply_bc(y + dy, Ny, bc.periodic_y());
        const int nz = apply_bc(z + dz, Nz, bc.periodic_z());

        if (nx < 0 || ny < 0 || nz < 0) return -1;
        return index_raw(nx, ny, nz);
    }

    void initialize_all(const Vec3& m_init) {
        switch(llg::init_type) {
            case llg::InitType::Uniform:
                for (auto& spin : spins) {
                    spin = m_init.normalized();
                }
                break;

            case llg::InitType::AFM_Checkerboard:
                for (int z = 0; z < Nz; ++z) {
                    for (int y = 0; y < Ny; ++y) {
                        for (int x = 0; x < Nx; ++x) {
                            int idx = index_raw(x, y, z);
                            if ((x + y + z) % 2 == 0) {
                                spins[static_cast<size_t>(idx)] = m_init.normalized();
                            } else {
                                spins[static_cast<size_t>(idx)] = (m_init * -1.0).normalized();
                            }
                        }
                    }
                }
                break;

            case llg::InitType::Random:
                for (auto& spin : spins) {
                    double theta = static_cast<double>(rand()) / RAND_MAX * 3.14159265358979323846;
                    double phi = static_cast<double>(rand()) / RAND_MAX * 2.0 * 3.14159265358979323846;
                    spin = Vec3(
                        std::sin(theta) * std::cos(phi),
                        std::sin(theta) * std::sin(phi),
                        std::cos(theta)
                    );
                }
                break;
        }
    }

    int index_raw(int x, int y, int z) const{
        return x + Nx * (y + Ny * z);
    }


    int Nx, Ny, Nz;
    Vec3i bc;

private:
    std::vector<Vec3> spins;

    static int apply_bc(int i, int N, bool periodic)
    {
        if (periodic) {
            i %= N;
            if (i < 0) i += N;
            return i;
        }
        return (i < 0 || i >= N) ? -1 : i;
    }
};


#endif // GRID_H