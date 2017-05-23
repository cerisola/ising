#ifndef THERMO_H
#define THERMO_H

typedef struct {
    double T;
    double J;
    double B;
    double transition_probabilities[2][5];
} Parameters;

int calculate_transition_probabilities(Parameters * parameters);
double get_transition_probability(int s, int sum_neighbours, const Parameters * parameters);

double magnetization(int * lattice, int n);
double energy(int * lattice, int n);

#endif //THERMO_H
