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
    if (argc < 11) {
        printf("usage: n Tmin Tmax J B npoints rounding nsep nsamples outdir (seed)\n");
        return 1;
    }
    unsigned int random_seed;
    if (argc == 12) {
        random_seed = atoi(argv[11]);
    } else {
        random_seed = (unsigned int)time(NULL);
    }

    int n = atoi(argv[1]);
    int * lattice = malloc(n * n * sizeof(int));

    double Tmin = atof(argv[2]);
    double Tmax = atof(argv[3]);
    Parameters parameters = { .T = Tmax, .J = atof(argv[4]), .B = atof(argv[5]) };

    int npoints = atoi(argv[6]);
    int rounding = atoi(argv[7]);
    int nsep = atoi(argv[8]) * n * n;
    int nsamples = atoi(argv[9]);
    int niter = (nsamples + 1) * nsep;

    pcg32_srand(random_seed);

    double prob = 0.5;
    fill_lattice(lattice, n, prob);

    ThermodynamicQuantities quantities;
    double * Mavg = malloc(npoints * sizeof(*Mavg));
    double * Eavg = malloc(npoints * sizeof(*Eavg));
    double * Aavg = malloc(npoints * sizeof(*Aavg));
    double * Mvar = malloc(npoints * sizeof(*Mvar));
    double * Evar = malloc(npoints * sizeof(*Evar));
    double * Avar = malloc(npoints * sizeof(*Avar));

    calculate_transition_probabilities(&parameters);
    set_thermodynamic_quantities(lattice, n, &parameters, &quantities);
    int ntherm = 10 * nsep;
    for (int i = 0; i < ntherm; i++) {
        metropolis(lattice, n, &parameters, &quantities);
    }

    double * Tvalues = create_linear_grid(Tmin, Tmax, npoints, rounding);
    double * Eval = malloc(nsamples * sizeof(*Eval));
    double * Mval = malloc(nsamples * sizeof(*Mval));
    double * Aval = malloc(nsamples * sizeof(*Aval));
    for (int j = 0; j < npoints; j++) {
        int i = npoints - 1 - j;
        parameters.T = Tvalues[i];
        calculate_transition_probabilities(&parameters);
        set_thermodynamic_quantities(lattice, n, &parameters, &quantities);

        Mavg[i] = 0;
        Eavg[i] = 0;
        Aavg[i] = 0;
        Mvar[i] = 0;
        Evar[i] = 0;
        Avar[i] = 0;
        double m = 0;
        double e = 0;
        double a = 0;
        int count = 0;
        for (int j = 0; j < niter; j++) {
            metropolis(lattice, n, &parameters, &quantities);
            m += quantities.M;
            e += quantities.E;
            a += fabs(quantities.M);
            count++;
            if (j > 0 && j % nsep == 0) {
                m /= count;
                e /= count;
                a /= count;
                update_online_mean_variance(m, j / nsep, Mavg + i, Mvar + i);
                update_online_mean_variance(e, j / nsep, Eavg + i, Evar + i);
                update_online_mean_variance(a, j / nsep, Aavg + i, Avar + i);
                Mval[j / nsep - 1] = m;
                Eval[j / nsep - 1] = e;
                Aval[j / nsep - 1] = a;
                m = 0;
                e = 0;
                a = 0;
                count = 0;
            }
        }
        Mvar[i] /= nsamples - 1;
        Evar[i] /= nsamples - 1;
        Avar[i] /= nsamples - 1;

        if (j > 0 && j % (npoints/10) == 0) {
            printf("Finished T_%d out of %d\n", j+1, npoints);
        }

        write_thermodynamic_quantities(argv[10], Mval, Eval, Aval, nsamples, n, &parameters, random_seed);
    }

    write_thermodynamic_quantities_temperature_sweep(argv[10],
            Tvalues, npoints, nsamples, Mavg, Eavg, Mvar, Evar, Aavg, Avar,
            nsep, n, &parameters, random_seed);

    free(Mavg);
    free(Eavg);
    free(Aavg);
    free(Mvar);
    free(Evar);
    free(Avar);
    free(Mval);
    free(Eval);
    free(Aval);
    free(Tvalues);
    free(lattice);

    return 0;
}
