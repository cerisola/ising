#include "metropolis.h"
#include <stdlib.h>
#include <math.h>
#include "random.h"

int metropolis(int * lattice, int n, const Parameters * parameters, ThermodynamicQuantities * quantities)
{
    int site = pick_site(lattice, n);
    flip(lattice, n, site, parameters, quantities);
    return 0;
}

int pick_site(int * lattice, int n)
{
    return pcg32_boundedrand(n * n);
}

int flip(int * lattice, int n, int site, const Parameters * parameters, ThermodynamicQuantities * quantities)
{
    int sum_neighbours = sum_neighbours_for_boundary_conditions(lattice, n, site, parameters->boundary_type);
    double p = get_transition_probability(lattice[site], sum_neighbours, parameters);
    if (p > 1 || pcg32_random() < p*PCG32_RAND_MAX) {
        lattice[site] = -lattice[site];
        update_thermodynamic_quantities(lattice[site], sum_neighbours, n, parameters, quantities);
    }
    return 0;
}
