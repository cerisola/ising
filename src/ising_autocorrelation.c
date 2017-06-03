#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>

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
     * outdir: output directory where data will be saved
     * seed: (optional) seed of the random number generator
    */
    if (argc < 6) {
        printf("usage: n T J B nsteps nsamples outdir (seed)\n");
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

    int nsteps = atoi(argv[5]);
    int nsamples = atoi(argv[6]);

    srand(random_seed);

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

    int niter = nsteps * (nsamples + 1);
    unsigned int accum = 0;
    for (int i = 0; i < niter; i++) {
        metropolis(lattice, n, &parameters, &quantities);

        if (i < nsteps) {
            int j = nsteps - 1 - i;
            Mval[j] = quantities.M;
            Eval[j] = quantities.E;
            Mcor[j] = 0;
            Ecor[j] = 0;
        } else {
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
        }

        if (i > 0 && i % (niter/100) == 0) {
            printf("Finished iter %d out of %d\n", i+1, niter);
        }
    }

    Mavg = Mavg / accum;
    Eavg = Eavg / accum;
    for (int i = 0; i < nsteps; i++) {
        Mcor[i] = Mcor[i] / accum - Mavg * Mavg;
        Ecor[i] = Ecor[i] / accum - Eavg * Eavg;
    }

    write_autocorrelation_values(argv[7], Mcor, Ecor, nsteps, nsamples, T, n, random_seed);

    free(Mval);
    free(Eval);
    free(Mcor);
    free(Ecor);
    free(lattice);

    return 0;
}
