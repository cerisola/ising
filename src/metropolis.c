#include "metropolis.h"
#include <stdlib.h>
#include <math.h>

int metropolis(int * lattice, int n, const Parameters * parameters, ThermodynamicQuantities * quantities)
{
    int site = pick_site(lattice, n);
    flip(lattice, n, site, parameters, quantities);
    return 0;
}

int pick_site(int * lattice, int n)
{
    return (int)(((double)rand())*n*n/RAND_MAX);
}

int flip(int * lattice, int n, int site, const Parameters * parameters, ThermodynamicQuantities * quantities)
{
    int i = site / n;
    int j = site % n;
    int sum_neighbours = lattice[((i+1) % n)*n + j] +
                         lattice[((i-1+n) % n)*n + j] +
                         lattice[i*n + ((j+1) % n)] +
                         lattice[i*n + ((j-1+n) % n)];
    double p = get_transition_probability(lattice[site], sum_neighbours, parameters);
    double q = ((double)rand())/RAND_MAX;
    if (q < p) {
        lattice[site] = -lattice[site];
        update_thermodynamic_quantities(lattice[site], sum_neighbours, n, parameters, quantities);
    }
    return 0;
}
