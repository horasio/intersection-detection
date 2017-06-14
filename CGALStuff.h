/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */

#ifndef RAY_CGAL_STUFF_H
#define RAY_CGAL_STUFF_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_3.h>
#include <CGAL/intersections.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
typedef K::Point_3 Point_3;
typedef K::Triangle_3 Triangle_3;
typedef K::Tetrahedron_3 Tetrahedron_3;

extern template class CGAL::Polyhedron_3<K>;
using Polyhedron_3 = CGAL::Polyhedron_3<K>;

#endif // RAY_CGAL_STUFF_H
