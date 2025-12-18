/**** Monte Carlo simulation of the 3D Ising model (anti/ferromagnet) with periodic boundary conditions ****/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define latticeSize 16
#define J 1.0

typedef struct {
    double x;
    double y;
    double z;
} vector3;

vector3 norm(vector3 v) {
    double length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    vector3 normalized = {v.x / length, v.y / length, v.z / length};
    return normalized;
}

vector3 lattice[latticeSize][latticeSize][latticeSize];

/**** Periodic boundary: small neighbor index ****/
int GetMinus(int index) {
    if (index == 0) return latticeSize - 1;
    return index - 1;
}

/**** Periodic boundary: high neighbor index ****/
int GetPlus(int index) {
    if (index == latticeSize - 1) return 0;
    return index + 1;
}

/**** Initialize lattice ****/
void fillLattice(void) {
    for (int i = 0; i < latticeSize; i++) {
        for (int j = 0; j < latticeSize; j++) {
            for (int k = 0; k < latticeSize; k++) {
                lattice[i][j][k] = (vector3){0.0, 0.0, 1.0};  // start ordered
            }
        }       
    }
}


int main(void) {
    return 0;
}