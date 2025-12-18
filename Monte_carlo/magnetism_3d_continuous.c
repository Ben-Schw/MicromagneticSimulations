/**** Monte Carlo simulation of the 3D Heisenberg model (anti/ferromagnet) with periodic boundary conditions ****/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define latticeSize 16
#define J 1.0
#define M_PI 3.14159265358979323846

/**** simple 3D vector structure ****/
typedef struct {
    double x;
    double y;
    double z;
} vector3;

/**** basic vector algebra ****/
double dot(vector3 a, vector3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vector3 norm(vector3 v) {
    double length = sqrt(dot(v, v));
    vector3 normalized = {v.x / length, v.y / length, v.z / length};
    return normalized;
}

/*** generate random unit vector ****/
vector3 random_unit_vector() {
    double theta = ((double)rand() / RAND_MAX) * 2.0 * M_PI;
    double phi = acos(1.0 - 2.0 * ((double)rand() / RAND_MAX));
    vector3 v = {sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi)};
    return v;
}

/*** 3D Lattice and magnetic field vector ***/
vector3 lattice[latticeSize][latticeSize][latticeSize]; // 3D lattice of spins

const vector3 B = {0.0, 0.0, 0.1}; // External magnetic field

/**** Periodic boundary: small neighbor index ****/
int getMinus(int index) {
    if (index == 0) return latticeSize - 1;
    return index - 1;
}

/**** Periodic boundary: high neighbor index ****/
int getPlus(int index) {
    if (index == latticeSize - 1) return 0;
    return index + 1;
}

/**** Magnetic Quantities****/
vector3 magentisation(void) {
    vector3 mag = {0.0, 0.0, 0.0};

    for (int i = 0; i < latticeSize; i++) {
        for (int j = 0; j < latticeSize; j++) {
            for (int k = 0; k < latticeSize; k++) {
                mag.x += (1.0 / (latticeSize * latticeSize *latticeSize)) * lattice[i][j][k].x;
                mag.y += (1.0 / (latticeSize * latticeSize *latticeSize)) * lattice[i][j][k].y;
                mag.z += (1.0 / (latticeSize * latticeSize *latticeSize)) * lattice[i][j][k].z;
            }
        }
    }

    return mag;
}

vector3 neelVector(void) {
    vector3 neel = {0.0, 0.0, 0.0};

    for (int i = 0; i < latticeSize; i++) {
        for (int j = 0; j < latticeSize; j++) {
            for (int k = 0; k < latticeSize; k++) {
                int stagger = ((i + j + k) % 2 == 0) ? 1 : -1;
                neel.x += (1.0 / (latticeSize * latticeSize *latticeSize)) * lattice[i][j][k].x * (double)stagger;
                neel.y += (1.0 / (latticeSize * latticeSize *latticeSize)) * lattice[i][j][k].y * (double)stagger;
                neel.z += (1.0 / (latticeSize * latticeSize *latticeSize)) * lattice[i][j][k].z * (double)stagger;
            }
        }
    }

    return neel;
}

/**** helper functions for simulation****/
vector3 nn_sum(int i, int j, int k){
    vector3 s = {0.0, 0.0, 0.0};
    vector3 a = {0.0, 0.0, 0.0};

    a = lattice[getMinus(i)][j][k];  s.x += a.x; s.y += a.y; s.z += a.z;
    a = lattice[getRight(i)][j][k]; s.x += a.x; s.y += a.y; s.z += a.z;
    a = lattice[i][getMinus(j)][k];  s.x += a.x; s.y += a.y; s.z += a.z;
    a = lattice[i][getRight(j)][k]; s.x += a.x; s.y += a.y; s.z += a.z;
    a = lattice[i][j][getMinus(k)];  s.x += a.x; s.y += a.y; s.z += a.z;
    a = lattice[i][j][getRight(k)]; s.x += a.x; s.y += a.y; s.z += a.z;

    return s;
}

/**** Hamiltonian difference ****/
double hamiltonian(int i, int j, int k, vector3 S_new) {
    vector3 S_old = lattice[i][j][k];
    vector3 sum_nn = nn_sum(i, j, k);
    vector3 diff = {S_new.x - S_old.x, S_new.y - S_old.y, S_new.z - S_old.z};

    double dH = -J * (dot(S_new, sum_nn) - dot(S_old, sum_nn)) - dot(B, diff);
    return dH;
}

/**** Metropolis algorithm  ****/
int metropolis(double temp) {
    for (int i = 0; i < latticeSize; i++) {
        for (int j = 0; j < latticeSize; j++) {
            for (int k = 0; k < latticeSize; k++) {

                vector3 S_new = random_unit_vector();
                double dH = hamiltonian(i, j, k, S_new);

                if (rand() < RAND_MAX * exp(-dH / temp)) {
                    lattice[i][j][k] = S_new;
                }
            }
        }
    }
    return 0;
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

/**** perform steps and write into output  ****/
void mc_steps(double temp, int limit, FILE *fp) {
    vector3 mag = {0.0, 0.0, 0.0};
    vector3 neel = {0.0, 0.0, 0.0};
    // thermalization and measurement
    for (int i = 0; i < limit; i++) {
        metropolis(temp);

        // after thermalization -> measure
        if (i >= 0.7 * limit) {
            vector3 m = magentisation();
            mag.x += m.x / (0.3 * limit);
            mag.y += m.y / (0.3 * limit);
            mag.z += m.z / (0.3 * limit);

            vector3 n = neelVector();
            neel.x += n.x / (0.3 * limit);
            neel.y += n.y / (0.3 * limit);
            neel.z += n.z / (0.3 * limit);
        }
    }
    // write output
    fprintf(fp, "%g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
            temp,
            B.x,
            B.y,
            B.z,
            mag.x,
            mag.y,
            mag.z,
            neel.x,
            neel.y,
            neel.z);
}

/**** main simulation ****/
int main(void) {
    srand((unsigned)time(NULL));

    // open output file
    FILE *f = fopen("3dLatticeSize16.csv", "w");
    if (!f) { perror("fopen"); return 1; }

    fprintf(f, "T, B_x, B_y, B_z, S_x, S_y, S_z, N_x, N_y, N_z\n");

    fillLattice();

    // simulation for several temperatures
    double temp;
    for (temp = 0.0; temp <= 4.0; temp += 0.02) {

        if (temp <= 2.0) {
            mc_steps(temp, 1000, f);
        } else if (temp < 2.4) {
            mc_steps(temp, 100000, f);
        } else {
            mc_steps(temp, 10000, f);
        }

    }

    fclose(f);
    return 0;
}