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
    int sum_neighbours = sum_neighbours_periodic_boundary_conditions(lattice, n, site);
    double p = get_transition_probability(lattice[site], sum_neighbours, parameters);
    if (p > 1 || rand() < p*RAND_MAX) {
        lattice[site] = -lattice[site];
        update_thermodynamic_quantities(lattice[site], sum_neighbours, n, parameters, quantities);
    }
    return 0;
}
