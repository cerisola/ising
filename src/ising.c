#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "metropolis.h"
#include "lattice.h"
#include "thermo.h"

int main(int argc, char ** argv)
{
    /* read input arguments; if none provided fallback to default values */
    if (argc < 6) {
        printf("usage: n T J B nsamples (seed)\n");
        return 1;
    }
    unsigned int random_seed;
    if (argc == 7) {
        random_seed = atoi(argv[6]);
    } else {
        random_seed = (unsigned int)time(NULL);
    }

    int n = atoi(argv[1]);
    int * lattice = malloc(n * n * sizeof(int));

    Parameters parameters = { .T = atof(argv[2]), .J = atof(argv[3]), .B = atof(argv[4]) };
    calculate_transition_probabilities(&parameters);

    int nsep = 20 * n*n;
    int nsamples = atoi(argv[5]);
    int niter = (nsamples + 1) * nsep;

    srand(random_seed);

    double prob = 0.5;
    fill_lattice(lattice, n, prob);

    ThermodynamicQuantities quantities;
    set_thermodynamic_quantities(lattice, n, &parameters, &quantities);
    double Mavg = 0;
    double Eavg = 0;
    for (int i = 0; i < niter; i++) {
        metropolis(lattice, n, &parameters, &quantities);
        if (i > 0 && i % nsep == 0) {
            Mavg += fabs(quantities.M);
            Eavg += quantities.E;
        }
    }
    Mavg /= nsamples;
    Eavg /= nsamples;

    print_lattice(lattice, n);
    printf("<M> = %f \t\t <E> = %f\n", Mavg, Eavg);

    free(lattice);

    return 0;
}
