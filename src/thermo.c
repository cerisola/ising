#include "thermo.h"

double magnetization(int *lattice, int n)
{
    double M = 0;
    for (int i = 0; i < n*n; i++) {
        M += lattice[i];
    }
    return M/(n*n);
}

int site_energy(int site, int *lattice, int n)
{
    int i = site / n;
    int j = site % n;
    int sum_neighbours = lattice[((i+1) % n)*n + j] +
                         lattice[((i-1+n) % n)*n + j] +
                         lattice[i*n + ((j+1) % n)] +
                         lattice[i*n + ((j-1+n) % n)];
    int siteE = - sum_neighbours * lattice[site];
    return siteE;
}

double energy(int *lattice, int n)
{
    double E = 0;
    for (int i = 0; i < n*n; i++) {
        E += site_energy(i, lattice, n);
    }
    E /= 2;
    return E/(n*n);
}
