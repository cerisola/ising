#ifndef THERMO_H
#define THERMO_H

int sum_neighbours_periodic_boundary_conditions(const int * lattice, int n, int site);

typedef struct {
    double T;
    double J;
    double B;
    double transition_probabilities[2][5];
} Parameters;

int calculate_transition_probabilities(Parameters * parameters);
double get_transition_probability(int s, int sum_neighbours, const Parameters * parameters);

typedef struct {
    double E;
    double M;
} ThermodynamicQuantities;

double magnetization(const int * lattice, int n);
double energy(const int * lattice, int n, const Parameters * parameters);
int set_thermodynamic_quantities(const int * lattice, int n, const Parameters * parameters, ThermodynamicQuantities * quantities);
int update_thermodynamic_quantities(const int new_site, int sum_neighbours, int n, const Parameters * parameters, ThermodynamicQuantities * quantities);

#endif //THERMO_H
