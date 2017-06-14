/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */

#ifndef RAY_K_MEANS_H
#define RAY_K_MEANS_H

#include "vec.h"
#include <vector>

void kMeans(std::vector<Vec3f> & pts, const int start, const int end, int & middle);

#endif // RAY_K_MEANS_H
