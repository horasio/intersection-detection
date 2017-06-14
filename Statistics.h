/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */
#ifndef RAY_STATISTICS_H
#define RAY_STATISTICS_H

#include <cstddef>
#include <cmath>
#include <limits>

template< typename T >
class Statistics {
  T accu, sqaccu, mini, submaxi, maxi;
  size_t count;

  public:

  void addSample(const T & s) {
    accu += s;
    sqaccu += s * s;
    if( mini > s ) mini = s;
    if( maxi < s ) {
         submaxi = maxi;
         maxi = s;
    } else if( (submaxi < s) && (s < maxi) ) {
         submaxi = s;
    }
    ++count;
  }
  Statistics() : accu(0), sqaccu(0), mini(std::numeric_limits<T>::max()),
     submaxi(std::numeric_limits<T>::lowest()),
     maxi(std::numeric_limits<T>::lowest()), count(0) {}
  size_t nbSamples() const { return count; }
  T min() const { return mini; }
  T max() const { return maxi; }
  T second_max() const { return submaxi; }
  T total() const { return accu; }
  T sqTotal() const { return sqaccu; }
  double mean() const { return static_cast<double>(accu)/static_cast<double>(count); }
  double dev() const {
    double temp = static_cast<double>(sqaccu)/static_cast<double>(count);
    double m = mean();
    return ::sqrt(temp - (m * m));
  }
};

#endif // RAY_STATISTICS_H
