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
     * outdir: output directory where data will be saved
     * seed: (optional) seed of the random number generator
    */
    if (argc < 8) {
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

    int nadiabatic_steps = 100;
    int ntherm = 100 * n * n;
    for (int i = 0; i < nadiabatic_steps; i++) {
        parameters.T = (nadiabatic_steps - i) * T;
        calculate_transition_probabilities(&parameters);
        set_thermodynamic_quantities(lattice, n, &parameters, &quantities);
        for (int i = 0; i < ntherm; i++) {
            metropolis(lattice, n, &parameters, &quantities);
        }
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
        double deltaM = fabs(quantities.M) - Mavg;
        double deltaE = quantities.E - Eavg;
        Mavg += deltaM / accum;
        Eavg += deltaE / accum;

        for (int j = 0; j < nsteps-1; j++) {
            Mval[j] = Mval[j+1];
            Eval[j] = Eval[j+1];
        }
        Mval[nsteps-1] = quantities.M;
        Eval[nsteps-1] = quantities.E;

        for (int j = 0; j < nsteps; j++) {
            double deltaM2 = fabs(Mval[j]) - Mavg;
            double deltaE2 = Eval[j] - Eavg;
            Mcor[j] += deltaM * deltaM2;
            Ecor[j] += deltaE * deltaE2;
        }

        if (i > 0 && i % (nsamples/100) == 0) {
            printf("Finished iter %d out of %d\n", i+1, nsamples);
        }
    }

    for (int i = 0; i < nsteps; i++) {
        Mcor[i] = Mcor[i] / accum;
        Ecor[i] = Ecor[i] / accum;
    }

    write_autocorrelation_values(argv[7], Mcor, Ecor, nsteps, nsamples, n, &parameters, random_seed);

    free(Mval);
    free(Eval);
    free(Mcor);
    free(Ecor);
    free(lattice);

    return 0;
}
