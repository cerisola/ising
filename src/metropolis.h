#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "thermo.h"

int metropolis(int * lattice, int n, const Parameters * parameters);
int pick_site(int * lattice, int n);
int flip(int * lattice, int n, int site, const Parameters * parameters);

#endif
