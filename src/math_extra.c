/* Copyright 2017 Federico Cerisola */
/* MIT License (see root directory) */
/* see header file for detailed documentation of each function */

#include <stdlib.h>
#include <math.h>

double * create_linear_grid(const double xmin, const double xmax,
                            const int npoints, const int round_digits)
{
    double slope = (xmax - xmin)/(npoints - 1);

    double * grid = malloc(npoints * sizeof(*grid));
    for (int n = 0; n < npoints; n++) {
        grid[n] = xmin + n * slope;
    }

    if (round_digits) {
        for (int n = 0; n < npoints; n++) {
            grid[n] = ((int)(grid[n]*((int)pow(10, round_digits))))/((double)pow(10, round_digits));
        }
    }

    return grid;
}

int update_online_mean_variance(double new_value, int n, double * mean, double * variance)
{
    double delta = new_value - (*mean);
    (*mean) += delta / n;

    double delta2 = new_value - (*mean);
    (*variance) += delta * delta2;

    return 0;
}
