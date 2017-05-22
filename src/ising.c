#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "metropolis.h"
#include "lattice.h"
#include "thermo.h"

int main(int argc, char ** argv)
{
    int n = 32;
    int * lattice = malloc(n * n * sizeof(int));
    double prob = 0.5;
    double T = 2.0;
    int nsamples = 2000;
    int nsep = 20 * n*n;
    int niter = (nsamples + 1) * nsep;
    srand(time(NULL));

    fill_lattice(lattice, n, prob);
    double Mavg = 0;
    double Eavg = 0;
    for (int i = 0; i < niter; i++) {
        metropolis(lattice, n, T);
        if (i > 0 && i % nsep == 0) {
            Mavg += fabs(magnetization(lattice, n));
            Eavg += energy(lattice, n);
        }
    }
    Mavg /= nsamples;
    Eavg /= nsamples;

    print_lattice(lattice, n);
    printf("<M> = %f \t\t <E> = %f\n", Mavg, Eavg);

    free(lattice);

    return 0;
}
