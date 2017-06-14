/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */

#ifndef _RAY_GJK_H_
#define _RAY_GJK_H_

#include "vec.h"
#include <iostream>

inline
Vec3f nearest_on_segment(const Vec3f & A, const Vec3f & B) {
  Vec3f AB( B - A );
  float u = A | AB;
  if( u >= 0.0f )
    return A;
  else {
    float l = AB.squaredLength();
    if( u <= -l )
      return B;
    else
      return A - ((u * AB) / l);
  }
}

struct GJK {
  Vec3f pts_[4];
  int n_;

  GJK() : n_(0) {}

  float delta_edge[6][2];
  float delta_triangle[4][3];
  float delta_tet[4];

  void clear() { n_ = 0; }

  bool closest_pt_on_tet(Vec3f & dir);

  void set(const Vec3f & p) {
    pts_[0] = p;
    n_ = 1;
  }

  bool add_and_nearest(const Vec3f & A, Vec3f & dir);
  bool add_and_nearest_for_distance(const Vec3f & A, Vec3f & dir);
  float add_and_distance(const Vec3f & A);
};

#endif // _RAY_GJK_H_
