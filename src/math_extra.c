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
