#include "metropolis.h"
#include <stdlib.h>
#include <math.h>

int metropolis(int * lattice, int n, float T)
{
    int site = pick_site(lattice, n);
    flip(lattice, n, T, site);
    return 0;
}

int pick_site(int * lattice, int n)
{
    return (int)(((double)rand())*n*n/RAND_MAX);
}

int flip(int * lattice, int n, float T, int site)
{
    int i = site / n;
    int j = site % n;
    int sum_neighbours = lattice[((i+1) % n)*n + j] +
                         lattice[((i-1+n) % n)*n + j] +
                         lattice[i*n + ((j+1) % n)] +
                         lattice[i*n + ((j-1+n) % n)];
    int deltaE = 2 * sum_neighbours * lattice[site];
    double p = exp(-deltaE/T);
    double q = ((double)rand())/RAND_MAX;
    if (q < p) {
        lattice[site] = -lattice[site];
    }
    return 0;
}
