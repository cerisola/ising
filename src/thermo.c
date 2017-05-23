#include "thermo.h"
#include <stdlib.h>
#include <math.h>

int calculate_transition_probabilities(Parameters * parameters)
{
    double deltaE_0 = -8 * parameters->J + 2 * parameters->B;
    double deltaE_1 = -4 * parameters->J + 2 * parameters->B;
    double deltaE_2 = 2 * parameters->B;
    double deltaE_3 = 4 * parameters->J + 2 * parameters->B;
    double deltaE_4 = 8 * parameters->J + 2 * parameters->B;
    parameters->transition_probabilities[0][0] = exp(+ deltaE_0 / parameters->T);
    parameters->transition_probabilities[0][1] = exp(+ deltaE_1 / parameters->T);
    parameters->transition_probabilities[0][2] = exp(+ deltaE_2 / parameters->T);
    parameters->transition_probabilities[0][3] = exp(+ deltaE_3 / parameters->T);
    parameters->transition_probabilities[0][4] = exp(+ deltaE_4 / parameters->T);
    parameters->transition_probabilities[1][0] = exp(- deltaE_0 / parameters->T);
    parameters->transition_probabilities[1][1] = exp(- deltaE_1 / parameters->T);
    parameters->transition_probabilities[1][2] = exp(- deltaE_2 / parameters->T);
    parameters->transition_probabilities[1][3] = exp(- deltaE_3 / parameters->T);
    parameters->transition_probabilities[1][4] = exp(- deltaE_4 / parameters->T);
    return 0;
}

double get_transition_probability(int s, int sum_neighbours, const Parameters * parameters)
{
    return parameters->transition_probabilities[(s + 1)/2][sum_neighbours/2 + 2];
}

double magnetization(int *lattice, int n)
{
    double M = 0;
    for (int i = 0; i < n*n; i++) {
        M += lattice[i];
    }
    return M/(n*n);
}

int site_interaction_energy(int site, int *lattice, int n, const Parameters * parameters)
{
    int i = site / n;
    int j = site % n;
    int sum_neighbours = lattice[((i+1) % n)*n + j] +
                         lattice[((i-1+n) % n)*n + j] +
                         lattice[i*n + ((j+1) % n)] +
                         lattice[i*n + ((j-1+n) % n)];
    int siteE = - parameters->J * sum_neighbours * lattice[site];
    return siteE;
}

double energy(int *lattice, int n, const Parameters * parameters)
{
    double E = 0;
    for (int i = 0; i < n*n; i++) {
        E += site_interaction_energy(i, lattice, n, parameters)/2 + parameters->B * lattice[i];
    }
    return E/(n*n);
}
