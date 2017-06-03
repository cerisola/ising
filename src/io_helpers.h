/*!
    @file io_helpers.h
    @brief This file contains helper functions for printing lattices.

    @author Federico Cerisola
    @copyright MIT License.
*/

#ifndef IO_HELPERS_H
#define IO_HELPERS_H

/*! Write lattice to file.

    @param path path to the folder were the lattice will be written. If the file
        exists it will be overwritten.
    @param lattice pointer to lattice to be written.
    @param L the linear size of the lattice.
*/
void write_lattice_to_file(const char * path, const int * lattice, int L);

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
    @param nsep separation of the samples (in units of Metropolis steps).
    @param L the linear size of the lattice.
    @param seed the random number generator seed used at the beginning of the
        probability sweep.
*/
void write_thermodynamic_quantities_temperature_sweep(const char * path,
        const double * temperature_grid, int grid_npoints, int nsamples,
        const double * Mavg, const double * Eavg, const double * Mvar,
        const double * Evar, int nsep, int L, unsigned int seed);

#endif /* IO_HELPERS_H */
