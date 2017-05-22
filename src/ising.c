#include <stdlib.h>
#include <stdio.h>
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
    int niter = 2000;
    srand(time(NULL));

    fill_lattice(lattice, n, prob);
    double Mavg = 0;
    double Eavg = 0;
    for (int i = 0; i < niter; i++) {
        metropolis(lattice, n, T);
        Mavg += magnetization(lattice, n);
        Eavg += energy(lattice, n);
    }
    Mavg /= niter;
    Eavg /= niter;

    print_lattice(lattice, n);
    printf("<M> = %f \t\t <E> = %f\n", Mavg, Eavg);

    free(lattice);

    return 0;
}
