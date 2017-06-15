/*!
    @file io_helpers.h
    @brief This file contains helper functions for printing lattices.

    @author Federico Cerisola
    @copyright MIT License.
*/

#ifndef IO_HELPERS_H
#define IO_HELPERS_H

#include "thermo.h"

/*! Write lattice to file.

    @param path path to the folder were the lattice will be written. If the file
        exists it will be overwritten.
    @param lattice pointer to lattice to be written.
    @param L the linear size of the lattice.
*/
void write_lattice_to_file(const char * path, const int * lattice, int L);

/*! Write thermodynamic quantities to file.

    @param path path to the folder where the data will be written. If the file
        exists it will be overwritten.
    @param Mval magnetization values.
    @param Eval energy values.
    @param Aux auxilliary additional magnitude to be saved. Pass NULL if no
        extra data should be saved.
    @param npoints number of data points. 
    @param L the linear size of the lattice.
    @param parameters struct with system parameters.
    @param seed the random number generator seed used at the beginning of the
        probability sweep.
*/
void write_thermodynamic_quantities(const char * path, const double * Mval,
        const double * Eval, const double * Aux, int npoints, int L,
        const Parameters * parameters, unsigned int seed);

/*! Write thermodynamic quantities statistics obtained in a temperature sweep to file.

    @param path path to the folder where the data will be written. If the file
        exists it will be overwritten.
    @param temperature_grid grid of the tested temperature values.
    @param grid_npoints the number of points in the temperature grid.
    @param nsampels the number of samples taken at each temperature.
    @param Mavg average values of the magnetization for each temperature.
    @param Eavg average values of the energy for each temperature.
    @param Mvar variance of the magnetization for each temperature.
    @param Evar variance of the energy of each temperature.
    @param Aavg average values of an additional auxilliarry magnitude. Pass
        NULL if none is desired to be saved.
    @param Avar variance of an additional auxilliary magnitude. Pass NULL if
        non is desired to be saved.
    @param nsep separation of the samples (in units of Metropolis steps).
    @param L the linear size of the lattice.
    @param parameters struct with system parameters.
    @param seed the random number generator seed used at the beginning of the
        probability sweep.
*/
void write_thermodynamic_quantities_temperature_sweep(const char * path,
        const double * temperature_grid, int grid_npoints, int nsamples,
        const double * Mavg, const double * Eavg, const double * Mvar,
        const double * Evar, const double * Aavg, const double * Avar,
        int nsep, int L, const Parameters * parameters, unsigned int seed);

/*! Write autocorrelation values to file.

    @param path path to the folder where the data will be written. If the file
        exists it will be overwritten.
    @param Mcor autocorrelation values of the magnetization.
    @param Ecor autocorrelation values of the energy.
    @param nsteps maximum number of steps used to calculate the autocorrelation.
    @param nsamples the number of samples taken at each step separation.
    @param L the linear size of the lattice.
    @param parameters struct with system parameters.
    @param seed the random number generator seed used at the beginning of the
        probability sweep.
*/
void write_autocorrelation_values(const char * path, const double * Mcor,
        const double * Ecor, int nsteps, int nsamples, int L,
        const Parameters * parameters, unsigned int seed);

#endif /* IO_HELPERS_H */
