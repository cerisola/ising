#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include "metropolis.h"
#include "lattice.h"
#include "thermo.h"
#include "math_extra.h"

int main(int argc, char ** argv)
{
    /* read input arguments; if none provided fallback to default values
     * n: lattice size
     * Tmin: minimum temperature value
     * Tmax: maximum temperature value
     * J: spin-spin coupling constant value
     * B: external magnetic field value
     * npoints: number of points in the temperature interval [Tmin, Tmax]
     * nsep: separation (in units of n^2) between different samples
     * nsamples: number of samples to average at each temperature
     * outdir: output directory where data will be saved
     * seed: (optional) seed of the random number generator
     * */
    if (argc < 6) {
        printf("usage: n Tmin Tmax J B npoints nsep nsamples outdir (seed)\n");
        return 1;
    }
    unsigned int random_seed;
    if (argc == 11) {
        random_seed = atoi(argv[10]);
    } else {
        random_seed = (unsigned int)time(NULL);
    }

    int n = atoi(argv[1]);
    int * lattice = malloc(n * n * sizeof(int));

    double Tmin = atof(argv[2]);
    double Tmax = atof(argv[3]);
    Parameters parameters = { .T = atof(argv[2]), .J = atof(argv[4]), .B = atof(argv[5]) };

    int npoints = atoi(argv[6]);
    int nsep = atoi(argv[7]) * n * n;
    int nsamples = atoi(argv[8]);
    int niter = (nsamples + 1) * nsep;

    srand(random_seed);

    double prob = 0.5;
    fill_lattice(lattice, n, prob);

    ThermodynamicQuantities quantities;
    double * Mavg = malloc(npoints * sizeof(*Mavg));
    double * Eavg = malloc(npoints * sizeof(*Eavg));
    double * Mvar = malloc(npoints * sizeof(*Mvar));
    double * Evar = malloc(npoints * sizeof(*Evar));

    double * Tvalues = create_linear_grid(Tmin, Tmax, npoints, 0);
    for (int i = 0; i < npoints; i++) {
        parameters.T = Tvalues[i];
        calculate_transition_probabilities(&parameters);
        set_thermodynamic_quantities(lattice, n, &parameters, &quantities);

        Mavg[i] = 0;
        Eavg[i] = 0;
        Mvar[i] = 0;
        Evar[i] = 0;
        for (int j = 0; j < niter; j++) {
            metropolis(lattice, n, &parameters, &quantities);
            if (j > 0 && j % nsep == 0) {
                Mavg[i] += fabs(quantities.M);
                Eavg[i] += quantities.E;
                Mvar[i] += quantities.M * quantities.M;
                Evar[i] += quantities.E * quantities.E;
            }
        }
        Mvar[i] = (Mvar[i] - (Mavg[i] * Mavg[i])/nsamples)/(nsamples - 1);
        Evar[i] = (Evar[i] - (Eavg[i] * Eavg[i])/nsamples)/(nsamples - 1);
        Mavg[i] /= nsamples;
        Eavg[i] /= nsamples;

        printf("Finished T_%d out of %d\n", i+1, npoints);
    }

    free(Mavg);
    free(Eavg);
    free(Mvar);
    free(Evar);
    free(Tvalues);
    free(lattice);

    return 0;
}
