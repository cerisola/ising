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


/*! Update the mean and variance for an online sampling process.

    This function updates the values of the sample mean and variance from an
    online sampling process.

    @param new_value the new value to be used to update statistics.
    @param n the number of the new sample.
    @param mean where the mean is stored. On function call it contains the old
        estimate of the mean and therefore on the first call should be initialized
        to 0. On function return in contains the updated value. The returned value is
        unnormalized and should be divided by the total number of observations at the end
        of the sampling process.
    @param variance where the variance is stored. On function call it contains the old
        estimate of the variance and therefore on the first call should be initialized
        to 0. On function return in contains the updated value. The returned value is
        unnormalized and should be divided by the total number of observations at the end
        of the sampling process.
*/
int update_online_mean_variance(double new_value, int n, double * mean, double * variance);

#endif /* MATH_EXTRA_H */
