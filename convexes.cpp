/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */

#include "convexes.h"
#include "kmeans.h"
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>

#include <queue>
#include <iostream>

using namespace std;
  
FullConvex::FullConvex(const Polyhedron_3 & ch) {
  static VH_to_id vhToId;
  vhToId.clear();
  int id (0);
  for( auto v = ch.vertices_begin(); v != ch.vertices_end(); ++v ) {
    const Point_3 & p = v->point();
    vertices_.emplace_back(p.x(), p.y(), p.z());
    vhToId.emplace(v, id++);
  }
  for( auto h = ch.halfedges_begin(); h != ch.halfedges_end(); ++h ) {
    int edgeV0 = vhToId[h->vertex()];
    int edgeV1 = vhToId[h->opposite()->vertex()];
    if( edgeV0 >= edgeV1 ) continue;
    int face1V = vhToId[h->next()->vertex()];
    int face0V = vhToId[h->opposite()->next()->vertex()];
    edges_.emplace_back(edgeV0, edgeV1, face0V, face1V);
  }
  for( auto f = ch.facets_begin(); f != ch.facets_end(); ++f ) {
    auto h = f->halfedge();
    const Point_3 & p0 = h->opposite()->vertex()->point(); Vec3f v0(p0.x(), p0.y(), p0.z());
    const Point_3 & p1 = h->vertex()->point();             Vec3f v1(p1.x(), p1.y(), p1.z());
    const Point_3 & p2 = h->next()->vertex()->point();     Vec3f v2(p2.x(), p2.y(), p2.z());
    Vec3f n = Vec3f::cross(v1-v0, v2-v0);
    n.normalize();
    facets_.emplace_back(n, n | v0);
  }
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

template< class NT >
NT radians_of_degrees(const NT d) {
	return d * NT(M_PI) / NT(180);
}

const unsigned char Frustum::edges[12][2] = {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
const unsigned char Frustum::faces[12][2] = {{0,4},{0,2},{0,5},{0,3},{4,1},{2,1},{5,1},{3,1},{4,3},{2,4},{5,2},{3,5}};

Frustum::Frustum() {
  RandomPointOnSphere rps;
  Matrix3f m;
  m.setRotation(rps.vec3f(), drand48() * 2.0 * M_PI);
  float near = 0.1f;
  float far = 100.0f;
  float halfFov = radians_of_degrees<float>(40.0f);
  float tanhf = ::tanf(halfFov);
  float nearx = near * tanhf;
  float ratio = 9.0f/16.0f;
  float neary = ratio * nearx;
  float farx = far * tanhf;
  float fary = ratio * farx;
  // Compute the center and radius of the largest inscribed sphere
  float f = fary; // f <- min(farx, fary)
  float l = ::sqrt(f*f + far*far);
  float r = fary * (l - f) / far;
  if( 2.0f * r + near > far) {
    r = 0.5f * (far - near);
  }
  radius_ = r;
  center_ = m.mulVec3(Vec3f(0.0f, 0.0f, 0.0f));
  // Now, compute the vertices
  vertices_[0] = m.mulVec3(Vec3f(-nearx, -neary, far-near-r));
  vertices_[1] = m.mulVec3(Vec3f(+nearx, -neary, far-near-r));
  vertices_[2] = m.mulVec3(Vec3f(+nearx, +neary, far-near-r));
  vertices_[3] = m.mulVec3(Vec3f(-nearx, +neary, far-near-r));
  vertices_[4] = m.mulVec3(Vec3f(-farx,  -fary,  -r));
  vertices_[5] = m.mulVec3(Vec3f(+farx,  -fary,  -r));
  vertices_[6] = m.mulVec3(Vec3f(+farx,  +fary,  -r));
  vertices_[7] = m.mulVec3(Vec3f(-farx,  +fary,  -r));
  lo_ = hi_ = vertices_[0];
  for( int i = 1; i < 8; ++i ) {
    lo_ = vecmin(lo_, vertices_[i]);
    hi_ = vecmax(hi_, vertices_[i]);
  }
  // compute planes equation
  makePlane(0, 0, 1, 3);// 0-1-2-3-0 : 0
  makePlane(1, 4, 7, 5);// 4-7-6-5-4 : 1
  makePlane(2, 1, 5, 2);// 1-5-6-2-1 : 2
  makePlane(3, 4, 0, 7);// 0-3-7-4-0 : 3
  makePlane(4, 4, 5, 0);// 4-5-1-0-4 : 4
  makePlane(5, 3, 2, 7);// 3-2-6-7-3 : 5
  for( int dir = 0; dir < 3; ++dir ) {
    // view along -dir axis
    // -X ==> plane YZ
    // -Y ==> plane ZX
    // -Z ==> plane XY
    int c0 = (dir+1) % 3;
    int c1 = (dir+2) % 3;
    int pos = 0;
    for( int e = 0; e < 12; ++e ) {
      float f0 = planes_[faces[e][0]][dir];
      float f1 = planes_[faces[e][1]][dir];
      if( f0 * f1 < 0.0f ) {
        Vec3f *v0(&vertices_[edges[e][0]]);
        Vec3f *v1(&vertices_[edges[e][1]]);
        if( f1 > 0.0f )
          std::swap(v0, v1);
        float nx = (*v1)[c1] - (*v0)[c1];
        float ny = (*v0)[c0] - (*v1)[c0];
        float n = ::sqrt(nx*nx+ny*ny);
        if( n < 1e-10 )
          continue;
        nx /= n;
        ny /= n;
        sil_[dir][pos].set(nx, ny, nx*(*v0)[c0] + ny*(*v0)[c1]);
        ++pos;
      }
    }
    sil_[dir][pos].set(0.0f, 0.0f, 0.0f);
  }
}

// An implementation of Greene's algorithm
// The silhouette 2D line equations are computed in the constructor
bool
Frustum::hitsAAB(const OBB & box) const {
  Vec3f boxLo = box.aaLo();
  Vec3f boxHi = box.aaHi();
  // First we test the bounding boxes:
  for( int dir = 0; dir < 3; ++dir ) {
    ++planeStatPerPair;
    if( hi_[dir] < boxLo[dir] ) return false;
    ++planeStatPerPair;
    if( lo_[dir] > boxHi[dir] ) return false;
  }
  // next we test the box against the side of the frustum:
  Vec3f m;
  for( int f = 0; f < 6; ++f ) {
    float mini;
    box.AAminimizeInDirection(planes_[f], mini, m);
    ++planeStatPerPair;
    if( maxes_[f] < mini ) return false;
  }
  // Finally we test in 2D orthogonally to each axis
  if( ! hitsAAB2D<1,2>(box) ) return false;
  if( ! hitsAAB2D<2,0>(box) ) return false;
  if( ! hitsAAB2D<0,1>(box) ) return false;
  return true;
}

void
Frustum::makePlane(int f, int a, int b, int c) {
  const Vec3f & va = vertices_[a];
  const Vec3f & vb = vertices_[b];
  const Vec3f & vc = vertices_[c];

  planes_[f] = Vec3f::cross(vb-va, vc-va).normalized();
  maxes_[f]  = planes_[f] | va;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

typedef CGAL::Min_sphere_of_points_d_traits_3<K, float> MS3;
typedef CGAL::Min_sphere_of_spheres_d<MS3> Min_sphere;

SortedVertices::SortedVertices(int n, float shift) {
  RandomPointOnSphere rps;
  Vec3f center, s(shift, 0.0f, 0.0f);
  for( int i = 0; i < n; ++i ) {
    Vec3f v = rps.vec3f()+s;
    center = center + v;
    vertices_.push_back(v);
  }
  center_ = (1.0f / n) * center;
  makeHierarchy();
}

void SortedVertices::makeHierarchy() {
  typedef std::queue<std::pair<int, int>> Queue;
  Queue queue;
  queue.emplace(0, nbVertices());
  //cerr << endl;
  while( ! queue.empty() ) {
    auto bounds = queue.front();
    //cerr << ' ' << nodes_.size() << ':' << (bounds.second - bounds.first);
    queue.pop();
    Min_sphere ms;
    for( int i = bounds.first; i < bounds.second; ++i ) {
      ms.insert(Point_3(vertices_[i].x(), vertices_[i].y(), vertices_[i].z()));
    }
    Node node;
    node.radius_ = ms.radius();
    //cerr << '(' << (bounds.second-bounds.first) << ")r=" << node.radius_ << "; ";
    auto coord = ms.center_cartesian_begin();
    float x = *coord++;
    float y = *coord++;
    float z = *coord++;
    node.center_.set(x,y,z);
    node.firstChild_ = -1;
    node.split_ = -1;
    if( bounds.second - bounds.first > 32 ) {
      kMeans(vertices_, bounds.first, bounds.second, node.split_);
      node.firstChild_ = nodes_.size() + queue.size() + 1;
      nodes_.push_back(node);
      queue.emplace(bounds.first, node.split_);
      queue.emplace(node.split_, bounds.second);
    }
    else
      nodes_.push_back(node);
  }
}

