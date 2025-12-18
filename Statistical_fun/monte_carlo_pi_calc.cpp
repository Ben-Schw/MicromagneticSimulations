/*
Project to calculate the value of Pi using the Monte Carlo method.
This is done by randomly generating points in a unit square and checking
how many fall within a quarter circle inscribed within that square.
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>

int main() {
    // Open output file
    FILE *f = fopen("monte_carlo_pi_output.csv", "w");
    if (!f) {
        perror("fopen");
        return 1;
    }

    fprintf(f, "# Estimated value of Pi using Monte Carlo method.\n");
    fprintf(f, "# Number of Samples, Pi Estimate\n");

    // Perform Monte Carlo simulation for different sample sizes
    for (int numSamples = 300; numSamples <= 900000; numSamples += 300) {
        // Counter for points inside the circle
        int insideCircle = 0;
        // Seed the random number generator
        std::srand(static_cast<unsigned int>(std::time(nullptr)));

        for (int i = 0; i < numSamples; ++i) {
            // Generate random point (x, y) in [0, 1) x [0, 1)
            double x = static_cast<double>(std::rand()) / RAND_MAX;
            double y = static_cast<double>(std::rand()) / RAND_MAX;

            // write generated point to file (optional)
            //fprintf(f, "%f %f\n", x, y);

            // Check if the point is inside the quarter circle
            if (x * x + y * y <= 1.0) {
                ++insideCircle;
            }
        }

        // Estimate the value of Pi
        double piEstimate = 4.0 * static_cast<double>(insideCircle) / numSamples;

        // Output the result
        fprintf(f, "%d, %f\n", numSamples, piEstimate);

    }
    fclose(f);

    return 0;
}