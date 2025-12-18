/**** Monte Carlo simulation of the 2D Ising model (anti/ferromagnet) with periodic boundary conditions****/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define latticeSize 16
#define bn 0.125 // parameter for finite-size scaling
#define gn 1.25
#define nuinvers 1
#define J 1.0
#define B 0.1

int lattice[latticeSize][latticeSize];

int fillLattice(void);
int Metropolis(double temp);
double Hamiltonian(int i, int j);
double Magnetisierung(void);
double Suszeptibilitaet(double temp, double wurzel, double magsquared);
double NeelVector(void);
void write(double limit, FILE *fp, double temp);
int GetLower(int index);
int GetHigher(int index);

/**** Initialize lattice ****/
int fillLattice(void) {
    for (int i = 0; i < latticeSize; i++) {
        for (int j = 0; j < latticeSize; j++) {
            lattice[i][j] = 1;  // start ordered
        }
    }
    return 0;
}

/**** Metropolis algorithm (one sweep in fixed i,j order) ****/
int Metropolis(double temp) {
    for (int i = 0; i < latticeSize; i++) {
        for (int j = 0; j < latticeSize; j++) {
            double dH = Hamiltonian(i, j);

            // Keep acceptance rule exactly as requested (your "point 4")
            if (rand() < RAND_MAX * exp(-dH / temp)) {
                lattice[i][j] = -1 * lattice[i][j];
            }
        }
    }
    return 0;
}

/**** Write measured observables ****/
void write(double limit, FILE *fp, double temp) {
    double mean_m2 = 0.0;      // <m^2> with m = |magnetization per spin|
    double mean_m = 0.0;       // <m>
    double mean_neel = 0.0;    // <|Neel|>

    for (double i = 0; i < limit; i += 1.0) {
        Metropolis(temp);

        // measurement after thermalization
        if (i >= 0.7 * limit) {
            double m = Magnetisierung();
            mean_m2   += (m * m) / (0.3 * limit);
            mean_m    += m / (0.3 * limit);
            mean_neel += NeelVector() / (0.3 * limit);
        }
    }

    double m_mean_sq = pow(mean_m, 2.0);        // <m>^2
    double neel_sq   = pow(mean_neel, 2.0);     // <|Neel|>^2

    double sus = Suszeptibilitaet(temp, mean_m2, m_mean_sq);
    double sqrt_m2 = sqrt(mean_m2);

    // CSV-like output (no trailing comma)
    fprintf(fp, "%g, %g, %g, %g, %g, %g, %g, %g, %g\n",
            -(temp - 2.269) / 2.269,
            B,
            m_mean_sq,
            sqrt_m2,
            sus,
            pow(latticeSize, bn),
            pow(latticeSize, gn),
            pow(latticeSize, nuinvers),
            neel_sq);
}

/**** Energy change for flipping spin (i,j): ΔE = 2 J s_ij * sum_nn s_nn + B * s_ij ****/
double Hamiltonian(int i, int j) {
    double ham = (double)2 * J * lattice[i][j] *
                 (lattice[GetLower(i)][j] + lattice[GetHigher(i)][j] +
                  lattice[i][GetLower(j)] + lattice[i][GetHigher(j)]) + (double) B * lattice[i][j];
    return ham;
}

/**** Periodic boundary: lower neighbor index ****/
int GetLower(int index) {
    if (index == 0) return latticeSize - 1;
    return index - 1;
}

/**** Periodic boundary: higher neighbor index ****/
int GetHigher(int index) {
    if (index == latticeSize - 1) return 0;
    return index + 1;
}

/**** Magnetization per spin (absolute value) ****/
double Magnetisierung(void) {
    double mag = 0.0;

    for (int i = 0; i < latticeSize; i++) {
        for (int j = 0; j < latticeSize; j++) {
            mag += (1.0 / (latticeSize * latticeSize)) * lattice[i][j];
        }
    }

    return fabs(mag);
}

/**** Susceptibility using χ = N/T * ( <m^2> - <m>^2 ), with m per spin ****/
double Suszeptibilitaet(double temp, double mean_m2, double mean_m_sq) {
    double sus = pow((double)latticeSize, 2.0) / temp * (mean_m2 - mean_m_sq);
    return sus;
}

/**** Neel (staggered) order parameter per spin (absolute value) ****/
double NeelVector(void) {
    double neel = 0.0;

    for (int i = 0; i < latticeSize; i++) {
        for (int j = 0; j < latticeSize; j++) {
            int stagger = ((i + j) % 2 == 0) ? 1 : -1;
            neel += (1.0 / (latticeSize * latticeSize)) * lattice[i][j] * (double)stagger;
        }
    }

    return fabs(neel);
}

int main(void) {
    srand((unsigned)time(NULL));

    FILE *f = fopen("QuelleLatticeSize128f.txt", "w");
    if (!f) { perror("fopen"); return 1; }

    fprintf(f, "T, B, <m>^2, sqrt(<m^2>), chi, L^b, L^g, L^(1/nu), <|Neel|>^2\n");

    fillLattice();

    double incr = 0.02;
    for (double temp = incr; temp < 4.0; temp += incr) {

        if (temp <= 2.0) {
            write(1000.0, f, temp);
        } else if (temp < 2.4) {
            write(100000.0, f, temp);
        } else {
            write(10000.0, f, temp);
        }
    }

    fclose(f);
    return 0;
}
