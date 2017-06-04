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
     * nsteps: maximum number of steps to measure correlation
     * nsamples: number of samples to calculate correlation averages
     * boundary: the type of boundary conditions
     * outdir: output directory where data will be saved
     * seed: (optional) seed of the random number generator
    */
    if (argc < 9) {
        printf("usage: n T J B nsteps nsamples boundary outdir (seed)\n");
        return 1;
    }
    unsigned int random_seed;
    if (argc == 10) {
        random_seed = atoi(argv[9]);
    } else {
        random_seed = (unsigned int)time(NULL);
    }

    int n = atoi(argv[1]);
    int * lattice = malloc(n * n * sizeof(int));

    double T = atof(argv[2]);
    Parameters parameters = { .T = T, .J = atof(argv[3]), .B = atof(argv[4]), .boundary_type = parse_boundary_type(argv[7]) };

    int nsteps = atoi(argv[5]);
    int nsamples = atoi(argv[6]);

    pcg32_srand(random_seed);

    double prob = 0.5;
    fill_lattice(lattice, n, prob);

    ThermodynamicQuantities quantities;
    double Mavg = 0;
    double Eavg = 0;
    double * Mval = malloc(nsteps * sizeof(*Mval));
    double * Eval = malloc(nsteps * sizeof(*Eval));
    double * Mcor = malloc(nsteps * sizeof(*Mcor));
    double * Ecor = malloc(nsteps * sizeof(*Ecor));

    calculate_transition_probabilities(&parameters);
    set_thermodynamic_quantities(lattice, n, &parameters, &quantities);
    int ntherm = 100 * n * n;
    for (int i = 0; i < ntherm; i++) {
        metropolis(lattice, n, &parameters, &quantities);
    }

    for (int i = 0; i < nsteps; i++) {
        metropolis(lattice, n, &parameters, &quantities);
        Mval[i] = quantities.M;
        Eval[i] = quantities.E;
        Mcor[i] = 0;
        Ecor[i] = 0;
    }

    unsigned int accum = 0;
    for (int i = 0; i < nsamples; i++) {
        metropolis(lattice, n, &parameters, &quantities);

        accum++;
        Mavg += quantities.M;
        Eavg += quantities.E;
        Mcor[0] += quantities.M * quantities.M;
        Ecor[0] += quantities.E * quantities.E;
        for (int j = 1; j < nsteps; j++) {
            Mcor[j] += Mval[j] * quantities.M;
            Ecor[j] += Eval[j] * quantities.E;
        }
        for (int j = 0; j < nsteps-1; j++) {
            Mval[j] = Mval[j+1];
            Eval[j] = Eval[j+1];
        }
        Mval[nsteps-1] = quantities.M;
        Eval[nsteps-1] = quantities.E;

        if (i > 0 && i % (nsamples/100) == 0) {
            printf("Finished iter %d out of %d\n", i+1, nsamples);
        }
    }

    Mavg = Mavg / accum;
    Eavg = Eavg / accum;
    for (int i = 0; i < nsteps; i++) {
        Mcor[i] = Mcor[i] / accum - Mavg * Mavg;
        Ecor[i] = Ecor[i] / accum - Eavg * Eavg;
    }

    write_autocorrelation_values(argv[8], Mcor, Ecor, nsteps, nsamples, T, n, random_seed);

    free(Mval);
    free(Eval);
    free(Mcor);
    free(Ecor);
    free(lattice);

    return 0;
}
