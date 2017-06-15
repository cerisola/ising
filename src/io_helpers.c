/* Copyright 2017 Federico Cerisola */
/* MIT License (see root directory) */
/* see header file for detailed documentation of each function */

#include "io_helpers.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <time.h>

void write_lattice_to_file(const char * path, const int * lattice, int L)
{
    time_t current_time = time(NULL);
    char * prefix = "lattice";
    size_t file_full_path_length = strlen(path) + strlen(prefix) + 160;
    char * file_full_path = malloc(file_full_path_length * sizeof(*file_full_path));
    sprintf(file_full_path, "%s/%s_%dx%d.csv", path, prefix, L, L);

    FILE * file_handler = fopen(file_full_path, "w");
    fprintf(file_handler, ";L:%d\n", L);
    fprintf(file_handler, ";date:%s", asctime(localtime(&current_time)));
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            fprintf(file_handler, "%d", lattice[i*L + j]);
            if (j < L - 1) {
                fputs(",", file_handler);
            }
        }
        fputs("\n", file_handler);
    }

    fclose(file_handler);
    free(file_full_path);
}

void write_thermodynamic_quantities(const char * path, const double * Mval,
        const double * Eval, const double * Aux, int npoints, int L,
        const Parameters * parameters, unsigned int seed)
{
    time_t current_time = time(NULL);
    char * prefix = "quantities";
    size_t file_full_path_length = strlen(path) + strlen(prefix) + 160;
    char * file_full_path = malloc(file_full_path_length * sizeof(*file_full_path));
    sprintf(file_full_path, "%s/%s_%dx%d_%.*e_%u.csv", path, prefix,
            L, L, DBL_DIG-1, parameters->T, seed);

    FILE * file_handler = fopen(file_full_path, "w");
    fprintf(file_handler, ";L:%d\n", L);
    fprintf(file_handler, ";seed:%u\n", seed);
    fprintf(file_handler, ";npoints:%d\n", npoints);
    fprintf(file_handler, ";T:%.*e\n", DBL_DIG-1, parameters->T);
    fprintf(file_handler, ";J:%.*e\n", DBL_DIG-1, parameters->J);
    fprintf(file_handler, ";B:%.*e\n", DBL_DIG-1, parameters->B);
    fprintf(file_handler, ";date:%s", asctime(localtime(&current_time)));

    for (int i = 0; i < npoints; i++) {
        if (Aux == NULL) {
            fprintf(file_handler, "%.*e,%.*e\n", DBL_DIG-1, Mval[i],
                    DBL_DIG-1, Eval[i]);
        } else {
            fprintf(file_handler, "%.*e,%.*e,%.*e\n", DBL_DIG-1, Mval[i],
                    DBL_DIG-1, Eval[i], DBL_DIG-1, Aux[i]);
        }
    }

    fclose(file_handler);
    free(file_full_path);
}

void write_thermodynamic_quantities_temperature_sweep(const char * path,
        const double * temperature_grid, int grid_npoints, int nsamples,
        const double * Mavg, const double * Eavg, const double * Mvar,
        const double * Evar, const double * Aavg, const double * Avar,
        int nsep, int L, const Parameters * parameters, unsigned int seed)
{
    time_t current_time = time(NULL);
    char * prefix = "temperature_sweep";
    size_t file_full_path_length = strlen(path) + strlen(prefix) + 160;
    char * file_full_path = malloc(file_full_path_length * sizeof(*file_full_path));
    sprintf(file_full_path, "%s/%s_%dx%d_%d_%u.csv", path, prefix,
            L, L, grid_npoints, seed);

    FILE * file_handler = fopen(file_full_path, "w");
    fprintf(file_handler, ";L:%d\n", L);
    fprintf(file_handler, ";seed:%u\n", seed);
    fprintf(file_handler, ";grid_npoints:%d\n", grid_npoints);
    fprintf(file_handler, ";nsamples:%d\n", nsamples);
    fprintf(file_handler, ";nsep:%d\n", nsep);
    fprintf(file_handler, ";J:%.*e\n", DBL_DIG-1, parameters->J);
    fprintf(file_handler, ";B:%.*e\n", DBL_DIG-1, parameters->B);
    fprintf(file_handler, ";date:%s", asctime(localtime(&current_time)));

    for (int i = 0; i < grid_npoints; i++) {
        if (Aavg == NULL || Avar == NULL) {
            fprintf(file_handler, "%.*e,%.*e,%.*e,%.*e,%.*e\n",
                    DBL_DIG-1, temperature_grid[i], DBL_DIG-1, Mavg[i],
                    DBL_DIG-1, Eavg[i], DBL_DIG-1, Mvar[i], DBL_DIG-1, Evar[i]);
        } else {
            fprintf(file_handler, "%.*e,%.*e,%.*e,%.*e,%.*e,%.*e,%.*e\n",
                    DBL_DIG-1, temperature_grid[i], DBL_DIG-1, Mavg[i],
                    DBL_DIG-1, Eavg[i], DBL_DIG-1, Mvar[i], DBL_DIG-1, Evar[i],
                    DBL_DIG-1, Aavg[i], DBL_DIG-1, Avar[i]);
        }
    }

    fclose(file_handler);
    free(file_full_path);
}

void write_autocorrelation_values(const char * path, const double * Mcor,
        const double * Ecor, int nsteps, int nsamples, int L,
        const Parameters * parameters, unsigned int seed)
{
    time_t current_time = time(NULL);
    char * prefix = "autocorrelation";
    size_t file_full_path_length = strlen(path) + strlen(prefix) + 160;
    char * file_full_path = malloc(file_full_path_length * sizeof(*file_full_path));
    sprintf(file_full_path, "%s/%s_%dx%d_%.*e_%u.csv", path, prefix,
            L, L, DBL_DIG-1, parameters->T, seed);

    FILE * file_handler = fopen(file_full_path, "w");
    fprintf(file_handler, ";L:%d\n", L);
    fprintf(file_handler, ";seed:%u\n", seed);
    fprintf(file_handler, ";nsteps:%d\n", nsteps);
    fprintf(file_handler, ";nsamples:%d\n", nsamples);
    fprintf(file_handler, ";T:%.*e\n", DBL_DIG-1, parameters->T);
    fprintf(file_handler, ";J:%.*e\n", DBL_DIG-1, parameters->J);
    fprintf(file_handler, ";B:%.*e\n", DBL_DIG-1, parameters->B);
    fprintf(file_handler, ";date:%s", asctime(localtime(&current_time)));

    for (int i = 0; i < nsteps; i++) {
        fprintf(file_handler, "%.*e,%.*e\n", DBL_DIG-1, Mcor[i], DBL_DIG-1, Ecor[i]);
    }

    fclose(file_handler);
    free(file_full_path);
}

