#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * Compute the distances between a center point and a list of atom positions.
 *
 * @param center A pointer to the coordinates of the center point (array of 3 doubles).
 * @param positions A pointer to the array of atom positions (array of 3 * n_positions doubles).
 * @param n_positions The number of atoms (length of the positions array divided by 3).
 * @param box A pointer to the dimensions of the simulation box (array of 3 doubles).
 * @param distances A pointer to the array where computed distances will be stored (array of n_positions doubles).
 */
void compute_distances(double *center, double *positions, int n_positions, double *box, double *distances) {
    for (int i = 0; i < n_positions; i++) {
        double dx = fabs(center[0] - positions[3 * i]);
        double dy = fabs(center[1] - positions[3 * i + 1]);
        double dz = fabs(center[2] - positions[3 * i + 2]);

        // Apply periodic boundary conditions
        if (dx > box[0] / 2.0) {
            dx = box[0] - dx;
        }
        if (dy > box[1] / 2.0) {
            dy = box[1] - dy;
        }
        if (dz > box[2] / 2.0) {
            dz = box[2] - dz;
        }

        distances[i] = sqrt(dx * dx + dy * dy + dz * dz);

        // Uncomment the following line for debugging each distance calculation
        // printf("Center: [%f, %f, %f], Position: [%f, %f, %f], Distance: %f\n", 
        //    center[0], center[1], center[2], positions[3 * i], positions[3 * i + 1], positions[3 * i + 2], distances[i]);
    }
}

/**
 * Compute the number of atoms within a cutoff radius for each center point.
 *
 * @param centers A pointer to the array of center coordinates (array of 3 * n_centers doubles).
 * @param n_centers The number of center points (length of the centers array divided by 3).
 * @param positions A pointer to the array of atom positions (array of 3 * n_positions doubles).
 * @param n_positions The number of atoms (length of the positions array divided by 3).
 * @param box A pointer to the dimensions of the simulation box (array of 3 doubles).
 * @param rcutoff The cutoff radius for counting atoms.
 * @param N_values A pointer to the array where the number of atoms within the cutoff radius for each center will be stored (array of n_centers doubles).
 */
void compute_for_z(double *centers, int n_centers, double *positions, int n_positions, double *box, double rcutoff, double *N_values) {
    // Allocate memory for distances
    double *distances = (double *)malloc(n_positions * sizeof(double));

    // Loop over each center
    for (int c = 0; c < n_centers; c++) {
        double *center = &centers[3 * c];
        // Compute distances from the current center to all atom positions
        compute_distances(center, positions, n_positions, box, distances);

        // Count the number of atoms within the cutoff radius
        int N_count = 0;
        for (int i = 0; i < n_positions; i++) {
            if (distances[i] < rcutoff) {
                N_count++;
            }
        }
        N_values[c] = N_count;

        // Uncomment the following line for debugging each center's N_count
        // printf("Center %d: N_count = %d\n", c, N_count);
    }
    // Free allocated memory
    free(distances);
}
