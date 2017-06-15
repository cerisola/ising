#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include "random.h"
#include "metropolis.h"
#include "lattice.h"
#include "thermo.h"
#include "math_extra.h"
#include "io_helpers.h"

int main(int argc, char ** argv)
{
    /* read input arguments; if none provided fallback to default values
     * n: lattice size
     * T: temperature value
     * J: spin-spin coupling constant value
     * B: external magnetic field value
     * npoints: number of points to measure
     * nsep: separation between samples (in units of n^2)
     * outdir: output directory where data will be saved
     * seed: (optional) seed of the random number generator
     * */
    if (argc < 8) {
        printf("usage: n T J B npoints nsep outdir (seed)\n");
        return 1;
    }
    unsigned int random_seed;
    if (argc == 9) {
        random_seed = atoi(argv[8]);
    } else {
        random_seed = (unsigned int)time(NULL);
    }

    int n = atoi(argv[1]);
    int * lattice = malloc(n * n * sizeof(int));

    double T = atof(argv[2]);
    Parameters parameters = { .T = T, .J = atof(argv[3]), .B = atof(argv[4]) };

    unsigned int npoints = atoi(argv[5]);
    unsigned int nsep = atoi(argv[6]);
    if (nsep == 0) {
        nsep = 1;
    }

    pcg32_srand(random_seed);

    double prob = 0.5;
    fill_lattice(lattice, n, prob);

    ThermodynamicQuantities quantities;
    unsigned int data_size = npoints / nsep;
    double * Mval = malloc(data_size * sizeof(*Mval));
    double * Eval = malloc(data_size * sizeof(*Eval));

    int nadiabatic_steps = 200;
    int ntherm = 100 * n * n;
    for (unsigned int i = 0; i < nadiabatic_steps; i++) {
        parameters.T = (nadiabatic_steps - i) * T;
        calculate_transition_probabilities(&parameters);
        set_thermodynamic_quantities(lattice, n, &parameters, &quantities);
        for (unsigned int i = 0; i < ntherm; i++) {
            metropolis(lattice, n, &parameters, &quantities);
        }
    }

    for (unsigned int i = 0; i < npoints; i++) {
        metropolis(lattice, n, &parameters, &quantities);

        if (i % nsep == 0) {
            Mval[i/nsep] = quantities.M;
            Eval[i/nsep] = quantities.E;
        }

        if (i > 0 && i % (npoints/10) == 0) {
            printf("Finished iter %d out of %d\n", i+1, npoints);
        }
    }

    write_thermodynamic_quantities(argv[7], Mval, Eval, data_size, n, &parameters, random_seed);

    free(Mval);
    free(Eval);
    free(lattice);

    return 0;
}
