#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <cmath>

#include "vector.h"
#include "parameters.h"
#include "magnetic_field.h"

// helper function:
static inline int idx(int x, int y, int z) {
    return x + llg::Nx * (y + llg::Ny * z);
}

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


int main() {
    // -- init grid --
    
}