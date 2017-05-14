#include "lattice.h"
#include <stdlib.h>
#include <stdio.h>

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
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            putchar(lattice[i*n + j] > 0 ? '+' : '-');
            if (j < n-1) {
                putchar(' ');
            } else {
                putchar('\n');
            }
        }
    }
    return 0;
}
