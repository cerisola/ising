#include "lattice.h"
#include <stdlib.h>

int fill_lattice(int * lattice, int n, float p)
{
    double q;
    for (int i = 0; i < n*n; i++) {
        q = ((double)rand())/RAND_MAX;
        if (q < p) {
            lattice[i] = 1;
        } else {
            lattice[i] = -1;
        }
    }
    return 0;
}

int print_lattice(int * lattice, int n)
{
  return 0;
}
