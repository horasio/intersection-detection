/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */

#include "gjk.h"
#include <limits>

//#define OUTPUT_DEBUG

static const int edge[6][2]        = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
//                                      0      1      2      3      4      5
static const int triangle[4][3]    = {{1,2,3}, {0,2,3}, {0,1,3}, {0,1,2}};

bool GJK::closest_pt_on_tet(Vec3f & dir) {
  float dots[4][4];
  int edge_co_neg[6] = {0,0,0,0,0,0};
  int v_co_neg[4] = {0,0,0,0};
  for( int i = 0; i < 4; ++i )
    for( int j = i; j < 4; ++j )
    {
      dots[i][j] = dots[j][i] = pts_[i] | pts_[j];
    }
  // Compute Deltas for each edge. 6*2 values to compute
  for( int eidx = 0; eidx < 6; ++eidx ) {
    int p0 = edge[eidx][0];
    int p1 = edge[eidx][1];
    if( (delta_edge[eidx][0] = dots[p1][p1] - dots[p1][p0]) <= 0.0f )
      ++v_co_neg[p1];
    if( (delta_edge[eidx][1] = dots[p0][p0] - dots[p0][p1]) <= 0.0f )
      ++v_co_neg[p0];
  }
  // Compute Deltas for each triangle. 4*3 values to compute
  for( int tidx = 0; tidx < 4; ++tidx ) {
    for( int v = 0; v < 3; ++v ) {
      int j  = triangle[tidx][v];
      int i0 = triangle[tidx][(v+1)%3];
      int i1 = triangle[tidx][(v+2)%3];
      if( i0 > i1 ) std::swap(i0, i1);
      int eidx = i0 + i1;
      if( i0 == 0 ) --eidx;
      int k = i0;
      delta_triangle[tidx][v] =
        delta_edge[eidx][0] * ( dots[i0][k] - dots[i0][j] ) +
        delta_edge[eidx][1] * ( dots[i1][k] - dots[i1][j] );
    }
  }
  // Compute Deltas for the tet. 4 values to compute
  for( int v = 0; v < 4; ++v ) {
    int k = triangle[v][0];
    delta_tet[v] =
      delta_triangle[v][0] * (dots[triangle[v][0]][k] - dots[triangle[v][0]][v] ) +
      delta_triangle[v][1] * (dots[triangle[v][1]][k] - dots[triangle[v][1]][v] ) +
      delta_triangle[v][2] * (dots[triangle[v][2]][k] - dots[triangle[v][2]][v] );
  }
  // --------------------------------------
  int nb_tet(0);
  for( int tidx = 0; tidx < 4; ++tidx ) {
    if( delta_tet[tidx] > 0.0f ) {
      ++nb_tet;
      continue;
    }
    // here delta_tet[tidx] is negative, so it is possible that the nearest
    // feature is the triangle tidx.
    int nb_tri(0);
    for( int vt = 0; vt < 3; ++vt ) {
      if( delta_triangle[tidx][vt] > 0.0f ) {
        ++nb_tri;
        continue;
      }
      // here delta_triangle[tidx][vt] is negative, so it is possible that the
      // nearest feature is the edge opposite to vertex vt in this triangle.
      int i0 = triangle[tidx][(vt+1)%3];
      int i1 = triangle[tidx][(vt+2)%3];
      int eidx = i0 + i1;
      if( i0 == 0 || i1 == 0 )
        --eidx;
      ++(edge_co_neg[eidx]);
    }
    if( 3 == nb_tri ) { // yes indeed, the triangle is the nearest feature
#ifdef OUTPUT_DEBUG
      std::cout << "NU "
        << delta_triangle[tidx][0] << "*(" << pts_[triangle[tidx][0]] << ") + "
        << delta_triangle[tidx][1] << "*(" << pts_[triangle[tidx][1]] << ") + "
        << delta_triangle[tidx][2] << "*(" << pts_[triangle[tidx][2]] << ") / "
        << (delta_triangle[tidx][0]+delta_triangle[tidx][1]+delta_triangle[tidx][2])
        << std::endl;
#endif
      dir =
        (delta_triangle[tidx][0] * pts_[triangle[tidx][0]] +
         delta_triangle[tidx][1] * pts_[triangle[tidx][1]] +
         delta_triangle[tidx][2] * pts_[triangle[tidx][2]])
        /
        (delta_triangle[tidx][0]+delta_triangle[tidx][1]+delta_triangle[tidx][2]);
      if( tidx < 3 )
        std::swap(pts_[tidx], pts_[3]);
      return false;
    }
  }
  if( nb_tet == 4 ) {
    return true; // origin inside tetrahedron
  }
  // ..............
  for( int eidx = 0; eidx < 6; ++eidx ) {
    if( (delta_edge[eidx][0] <= 0.0f) || (delta_edge[eidx][1] <= 0.0f) || (edge_co_neg[eidx] < 2) )
      continue;
    Vec3f p0 = pts_[edge[eidx][0]];
    Vec3f p1 = pts_[edge[eidx][1]];
#ifdef OUTPUT_DEBUG
      std::cout << "NU "
        << delta_edge[eidx][0] << "*(" << p0 << ") + "
        << delta_edge[eidx][1] << "*(" << p1 << ") / "
        << (delta_edge[eidx][0]+delta_edge[eidx][1])
        << std::endl;
#endif
    dir = (delta_edge[eidx][0] * p0 + delta_edge[eidx][1] * p1)
      /
      (delta_edge[eidx][0]+delta_edge[eidx][1]);
    pts_[0] = p0;
    pts_[1] = p1;
    n_ = 2;
    return false;
  }
  for( int i = 0; i < 4; ++i ) {
    if( v_co_neg[i] < 3 ) continue;
#ifdef OUTPUT_DEBUG
      std::cout << "NU " << pts_[i] << std::endl;
#endif
    dir = pts_[i];
    if( i > 0 )
      pts_[0] = pts_[i];
    n_ = 1;
    return false;
  }
  // THE BACKUP PROCEDURE
  int itri = -1;
  int iedge = -1;
  int ivertex = -1;
  float dirn2 = std::numeric_limits<float>::max();
  for( int tidx = 0; tidx < 4; ++tidx ) {
    if( (delta_triangle[tidx][0] <= 0.0f) || (delta_triangle[tidx][1] <= 0.0f) || (delta_triangle[tidx][2] <= 0.0f) )
      continue;
    Vec3f candidate =
      (delta_triangle[tidx][0] * pts_[triangle[tidx][0]] +
       delta_triangle[tidx][1] * pts_[triangle[tidx][1]] +
       delta_triangle[tidx][2] * pts_[triangle[tidx][2]])
      /
      (delta_triangle[tidx][0]+delta_triangle[tidx][1]+delta_triangle[tidx][2]);
    float n2(candidate.norm2());
    if( dirn2 > n2 ) {
      dir = candidate;
#ifdef OUTPUT_DEBUG
       std::cout << "NU IN BACKUP" << std::endl;
#endif
      dirn2 = n2;
      itri = tidx;
    }
  }
  for( int eidx = 0; eidx < 6; ++eidx ) {
    if( (delta_edge[eidx][0] <= 0.0f) || (delta_edge[eidx][1] <= 0.0f) )
      continue;
    Vec3f candidate =
      (delta_edge[eidx][0] * pts_[edge[eidx][0]]+ delta_edge[eidx][1] * pts_[edge[eidx][1]])
      /
      (delta_edge[eidx][0]+delta_edge[eidx][1]);
    float n2(candidate.norm2());
    if( dirn2 > n2 ) {
      dir = candidate;
#ifdef OUTPUT_DEBUG
       std::cout << "NU IN BACKUP" << std::endl;
#endif
      dirn2 = n2;
      iedge = eidx;
    }
  }
  for( int v = 0; v < 4; ++v ) {
    if( dirn2 > dots[v][v] ) {
      dirn2 = dots[v][v];
      ivertex = v;
    }
  }
  if( ivertex >= 0 ) {
    dir = pts_[0] = pts_[ivertex];
#ifdef OUTPUT_DEBUG
       std::cout << "NU IN BACKUP" << std::endl;
#endif
    n_ = 1;
  } else if( iedge >= 0 ) {
    pts_[0] = pts_[edge[iedge][0]];
    pts_[1] = pts_[edge[iedge][1]];
    n_ = 2;
  } else if( itri >= 0 ) {
    if( itri < 3 )
      std::swap(pts_[itri], pts_[3]);
  } else
    std::cerr << "BAAAAAAAAAAAAAAAAAD";
  return false;
}

bool GJK::add_and_nearest(const Vec3f & A, Vec3f & dir) {
  // return true if the updated simplex contains the origin
  switch( n_ ) {
    case 0: {
              pts_[0] = A;
              n_ = 1;
              dir = A;
#ifdef OUTPUT_DEBUG
       std::cout << "NU FROM NOTHING " << dir << std::endl;
#endif
              return false;
              break;
            }
    case 1: {
              pts_[1] = A;
              n_ = 2;
              Vec3f AB = pts_[0] - A;
              dir = A - (((A | AB) / AB.squaredLength()) * AB);
              return false;
              break;
            }
    case 2: {
              Vec3f AO = A; AO.negate();
              Vec3f AB = pts_[0] - A;
              Vec3f AC = pts_[1] - A;
              Vec3f N = Vec3f::cross(AB, AC); // normal
              Vec3f outAB = Vec3f::cross(AB, N);
              if( (outAB | AO) >= 0.0f ) {
                pts_[1] = A;
                dir = Vec3f::cross(Vec3f::cross(AB, AO), AB);
#ifdef OUTPUT_DEBUG
       std::cout << "NU FROM SEG, out of AB " << dir << std::endl;
#endif
                return false;
              }
              Vec3f outAC = Vec3f::cross(N, AC);
              if( (outAC | AO) >= 0.0f ) {
                pts_[0] = A;
                dir = Vec3f::cross(Vec3f::cross(AC, AO), AC);
                return false;
              }
              pts_[2] = A;
              n_ = 3;
              dir = N;
              // We want the three points in pts_ to be ClockWise as seen from the origin
              if( (A | N) < 0.0f ) {
                std::swap(pts_[0], pts_[1]);
                dir.negate();
              }
              return false;
              break;
            }
    case 3: {
              pts_[3] = A;
              return closest_pt_on_tet(dir);
              break;
            }
    default: std::cerr << "BIG MISBEHAVIOR in GJK::nearest"; break;
  }
  return false;
}

bool GJK::add_and_nearest_for_distance(const Vec3f & A, Vec3f & dir) {
  // return true if the updated simplex contains the origin
  switch( n_ ) {
    case 0: {
              pts_[0] = A;
              n_ = 1;
              dir = A;
#ifdef OUTPUT_DEBUG
       std::cout << "NUFD FROM NOTHING " << dir << std::endl;
#endif
              return false;
              break;
            }
    case 1: {
              pts_[1] = A;
              n_ = 2;
              dir = nearest_on_segment(A, pts_[0]);
#ifdef OUTPUT_DEBUG
       std::cout << "NUFD FROM POINT " << dir << std::endl;
#endif
              return false;
              break;
            }
    case 2: {
              Vec3f AO = A; AO.negate();
              Vec3f AB = pts_[0] - A;
              Vec3f AC = pts_[1] - A;
              Vec3f N = Vec3f::cross(AB, AC); // normal
              Vec3f outAB = Vec3f::cross(AB, N);
              if( (outAB | AO) >= 0.0f ) {
                pts_[1] = A;
                dir = nearest_on_segment(A, pts_[0]);
                //dir = Vec3f::cross(Vec3f::cross(AB, AO), AB);
#ifdef OUTPUT_DEBUG
       std::cout << "NUFD FROM SEG, out of AB " << dir << std::endl;
#endif
                return false;
              }
              Vec3f outAC = Vec3f::cross(N, AC);
              if( (outAC | AO) >= 0.0f ) {
                pts_[0] = A;
                dir = nearest_on_segment(A, pts_[1]);
                //dir = Vec3f::cross(Vec3f::cross(AC, AO), AC);
#ifdef OUTPUT_DEBUG
       std::cout << "NUFD FROM SEG, out of AC " << dir << std::endl;
#endif
                return false;
              }
              pts_[2] = A;
              n_ = 3;
              dir = N;
              // We want the three points in pts_ to be ClockWise as seen from the origin
              if( (A | N) < 0.0f ) {
                std::swap(pts_[0], pts_[1]);
                dir.negate();
              }
              dir = (A | dir) * dir / dir.norm2();
#ifdef OUTPUT_DEBUG
       std::cout << "NUFD FROM SEG, NORMAL " << dir << std::endl;
#endif
              return false;
              break;
            }
    case 3: {
              pts_[3] = A;
              return closest_pt_on_tet(dir);
              break;
            }
    default: std::cerr << "BIG MISBEHAVIOR in GJK::nearest"; break;
  }
  return false;
}

float GJK::add_and_distance(const Vec3f & A) {
  Vec3f dir;
  add_and_nearest_for_distance(A, dir);
  dir.normalize();
  float m(dir | pts_[0]);
  for( int i = 1; i < n_; ++i )
    m = std::max(m, dir | pts_[i]);
  return m;
}

