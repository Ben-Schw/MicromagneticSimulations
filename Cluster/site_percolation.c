// site percolation simulation in C using C17 standard and Hoshen-Kopelman algorithm
// to identify clusters in a 2D lattice with periodic boundary conditions and nearest-neighbor connectivity.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define SIZE 10  // Lattice size
#define MAX_LABELS (SIZE * SIZE)

static int lattice[SIZE][SIZE];
static int parents[MAX_LABELS + 1];
static int cluster_sizes[MAX_LABELS + 1];

// Find with path compression
static int find_root(int x) {
    while (parents[x] != x) {
        parents[x] = parents[parents[x]];
        x = parents[x];
    }
    return x;
}

// Union by size
static int union_clusters(int x, int y) {
    int rootX = find_root(x);
    int rootY = find_root(y);

    if (rootX == rootY) return rootX;

    if (cluster_sizes[rootX] < cluster_sizes[rootY]) {
        int tmp = rootX;
        rootX = rootY;
        rootY = tmp;
    }

    parents[rootY] = rootX;
    cluster_sizes[rootX] += cluster_sizes[rootY];
    return rootX;
}

static void print_lattice(void) {
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            printf("%3d ", lattice[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Initialize lattice and union-find arrays
static void initialize(double p) {
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            lattice[i][j] = ((double)rand() / (double)RAND_MAX < p) ? 1 : 0;
        }
    }

    for (int i = 0; i <= MAX_LABELS; i++) {
        parents[i] = i;
        cluster_sizes[i] = 0;
    }
}

// Hoshen-Kopelman labeling (4-neighbor, without PBC unions here)
static int hoshen_kopelman(void) {
    int next_label = 1;

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            if (lattice[i][j] == 0) continue; // empty

            int left = (j > 0) ? lattice[i][j - 1] : 0;
            int up   = (i > 0) ? lattice[i - 1][j] : 0;

            if (left == 0 && up == 0) {
                // new cluster label
                int lbl = next_label++;
                lattice[i][j] = lbl;
                parents[lbl] = lbl;
                cluster_sizes[lbl] = 1; // count this site exactly once
            } else if (left != 0 && up == 0) {
                int rootL = find_root(left);
                lattice[i][j] = rootL;
                cluster_sizes[rootL] += 1; // count this site
            } else if (left == 0 && up != 0) {
                int rootU = find_root(up);
                lattice[i][j] = rootU;
                cluster_sizes[rootU] += 1; // count this site
            } else {
                // both neighbors occupied
                int root = union_clusters(left, up);
                lattice[i][j] = root;
                cluster_sizes[find_root(root)] += 1; // count this site once (after union)
            }
        }
    }

    return next_label - 1;
}

static void relabel_to_roots(void) {
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            if (lattice[i][j] != 0) lattice[i][j] = find_root(lattice[i][j]);
        }
    }
}

// Apply periodic boundary conditions: union left-right and top-bottom edges
static void apply_pbc(void) {
    for (int i = 0; i < SIZE; i++) {
        int a = lattice[i][0];
        int b = lattice[i][SIZE - 1];
        if (a != 0 && b != 0) union_clusters(a, b);
    }
    for (int j = 0; j < SIZE; j++) {
        int a = lattice[0][j];
        int b = lattice[SIZE - 1][j];
        if (a != 0 && b != 0) union_clusters(a, b);
    }
}

// Get maximum cluster size among roots
static int max_cluster_size(int max_label) {
    int max = 0;
    for (int lbl = 1; lbl <= max_label; lbl++) {
        if (parents[lbl] == lbl) {
            if (cluster_sizes[lbl] > max) max = cluster_sizes[lbl];
        }
    }
    return max;
}

int main(void) {
    srand((unsigned)time(NULL));

    FILE *out = fopen("ratio_vs_p.txt", "w");
    if (!out) { perror("fopen"); return 1; }

    const int runs = 10;

    for (double p = 0.0; p <= 1.0 + 1e-9; p += 0.1) {
        double average_ratio = 0.0;

        for (int r = 0; r < runs; r++) {
            initialize(p);

            if (fabs(p - 0.5) < 1e-9 && r == 0) {
                printf("p = 0.5, initial lattice (1=occupied,0=empty):\n");
                print_lattice();
            }

            int max_label = hoshen_kopelman();
            apply_pbc();
            relabel_to_roots();

            if (fabs(p - 0.5) < 1e-9 && r == 0) {
                printf("p = 0.5, labeled lattice (roots):\n");
                print_lattice();
            }

            int smax = max_cluster_size(max_label);
            average_ratio += (double)smax / (double)(SIZE * SIZE) / (double)runs;
        }

        fprintf(out, "%.2f %.6f\n", p, average_ratio);
        printf("p=%.1f  <Smax/N>=%.4f\n", p, average_ratio);
    }

    fclose(out);
    return 0;
}
