/*!
    @file math_extra.h
    @brief This file contains generic math functions not included in the
    standard library.

    @author Federico Cerisola
    @copyright MIT License.
*/

#ifndef MATH_EXTRA_H
#define MATH_EXTRA_H

/*! Create an equally spaced grid in a given interval.

    This function creates a linear discrete grid of \f$N\f$ points in the
    interval \f$[x_{min}, x_{max}]\f$.
    The distribution of points is given by the formula
    \f[
        x_n = x_min + n \frac{\left(x_max - x_min)}{N - 1}, \quad n = 0,\cdots,N-1
    \f]

    @param xmin initial value of the interval.
    @param xmax final value of the interval.
    @param npoints the number of points in the grid.
    @param round_digits to how many digits should the grid points be rounded.
        If round_digits = 0, then no rounding takes place.

    @returns pointer to the grid.

    @warning The returned grid must later be freed by the user.
*/
double * create_linear_grid(const double xmin, const double xmax,
                            const int npoints, const int round_digits);

#endif /* MATH_EXTRA_H */
