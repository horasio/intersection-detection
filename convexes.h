/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */

#ifndef RAY_CONVEXES_H
#define RAY_CONVEXES_H

#include "vec.h"
#include <algorithm>
#include <vector>
#include <stack>
#include <map>
#include "obb.h"
#include "Matrix3.h"
#include "CGALStuff.h"

extern unsigned long planeStatPerPair;

// ---------------------------------------------------------------

struct RandomPointOnSphere {
  Point_3 point3() const {
    float z = 2.0f * drand48() - 1.0f;
    float angle = 2.0f * M_PI * drand48();
    float rad = ::sqrtf(1.0f - z*z);
    return Point_3(rad * cosf(angle), rad * sinf(angle), z);
  }
  Vec3f vec3f() const {
    float z = 2.0f * drand48() - 1.0f;
    float angle = 2.0f * M_PI * drand48();
    float rad = ::sqrtf(1.0f - z*z);
    return Vec3f(rad * cosf(angle), rad * sinf(angle), z);
  }
  RandomPointOnSphere & operator++(int) {
    return (*this);
  }
  RandomPointOnSphere & operator++() {
    return (*this);
  }
};

// ---------------------------------------------------------------

struct alignas(16) Sphere {
  Vec3f center_;
  float radius_;

  Sphere(const Vec3f & c, const float r) : center_(c), radius_(r) {}

  const Vec3f & center() const { return center_; }

  Vec3f vertex(int) const { return Vec3f(center_.x()+radius_, center_.y(), center_.z()); }

  void maximizeInDirection(const Vec3f & dir, float & maxi, Vec3f & maximizer) const {
    float l = dir.length();
    maxi = (dir | center_) + radius_ * l;
    maximizer = center_ + ((radius_ / l) * dir);
  }
  void minimizeInDirection(const Vec3f & dir, float & mini, Vec3f & minimizer) const {
    float l = dir.length();
    mini = (dir | center_) - radius_ * l;
    minimizer = center_ - ((radius_ / l) * dir);
  }
};

// ---------------------------------------------------------------

struct alignas(16) Tet {
  Vec3f vertices_[4]; // positive orientation
  Vec3f normals_[4]; // 123, 032, 013, 021
  Vec3f center_;

  const Vec3f & vertex(const int i) const { return vertices_[i]; }

  inline
  void computeNormalsAndCenter() {
    normals_[0]=Vec3f::cross(vertices_[2]-vertices_[1], vertices_[3]-vertices_[1]);
    normals_[1]=Vec3f::cross(vertices_[3]-vertices_[0], vertices_[2]-vertices_[0]);
    normals_[2]=Vec3f::cross(vertices_[1]-vertices_[0], vertices_[3]-vertices_[0]);
    normals_[3]=Vec3f::cross(vertices_[2]-vertices_[0], vertices_[1]-vertices_[0]);
    center_.set(0.0f, 0.0f, 0.0f);
    for( int i = 0; i < 4; ++i ) {
      normals_[i].normalize();
      center_ = center_ + vertices_[i];
    }
    center_ = 0.25f * center_;
  }

  const Vec3f & center() const { return center_; }

  inline
  float minimizePlane(const Vec4f & plane) const {
    return
    std::min(
        std::min( plane(vertices_[0]), plane(vertices_[1])),
        std::min( plane(vertices_[2]), plane(vertices_[3])));
  }

  float minAlongDirection(const Vec3f & dir) const {
    float mm(vertices_[0] | dir);
    for( int i = 1; i < 4; ++i ) {
      float v = vertices_[i] | dir;
      if( v < mm )
        mm = v;
    }
    return mm;
  }

  float maxAlongDirection(const Vec3f & dir) const {
    float mm(vertices_[0] | dir);
    for( int i = 1; i < 4; ++i ) {
      float v = vertices_[i] | dir;
      if( v > mm )
        mm = v;
    }
    return mm;
  }

  void shift(float shift) {
    for( int i = 0; i < 4; ++i ) {
      vertices_[i].x() += shift;
    }
    center_.x() += shift;
  }

  bool containsZero() const {
    for( int i = 0; i < 4; ++i ) {
      if( (normals_[i] | vertices_[(i+1) % 4]) < 0.0f )
        return false; // tet does NOT contain the origin
    }
    return true;
  }

  bool faceSeparate(const Tet & t) const {
    for( int i = 0; i < 4; ) {
      Vec4f plane(normals_[i]);
      plane.w() = - ( plane | vertices_[(++i) & 0x3] );
      ++planeStatPerPair;
      if( t.minimizePlane(plane) >= 0.0f )
        return true; // tets do NOT intersect
    }
    return false;
  }

  inline
  bool differenceCoversZeroInDir(const Tet & B, int & vA, int & vB, const Vec3f & dir) const {
    // difference above is: A - B
    ++planeStatPerPair;
    float maxOverA = vertices_[0] | dir;
    float minOverB = B.vertices_[0] | dir;
    vA = vB = 0;
    for( int i = 1; i < 4; ++i ) {
      float tempA = vertices_[i] | dir;
      float tempB = B.vertices_[i] | dir;
      if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
      if( tempB < minOverB ) { minOverB = tempB; vB = i; }
    }
    return maxOverA >= minOverB;
  }

  inline
  void differenceInDir(const Tet & B, int & vA, int & vB, float & maxOverA, float & minOverB, const Vec3f & dir) const {
    // difference above is: A - B
    ++planeStatPerPair;
    maxOverA = vertices_[0] | dir;
    minOverB = B.vertices_[0] | dir;
    vA = vB = 0;
    for( int i = 1; i < 4; ++i ) {
      float tempA = vertices_[i] | dir;
      float tempB = B.vertices_[i] | dir;
      if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
      if( tempB < minOverB ) { minOverB = tempB; vB = i; }
    }
  }

  bool differenceCoversZeroInDirLazy(const Tet & B, int & vA, int & vB, const Vec3f & dir) const {
    // difference above is: A - B
    ++planeStatPerPair;
    float maxOverA = vertices_[0] | dir;
    float minOverB = B.vertices_[0] | dir;
    vA = vB = 0;
    if( maxOverA >= minOverB ) return true;
    for( int i = 1; i < 4; ++i ) {
      float tempA = vertices_[i] | dir;
      float tempB = B.vertices_[i] | dir;
      if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
      if( tempB < minOverB ) { minOverB = tempB; vB = i; }
      if( maxOverA >= minOverB ) return true;
    }
    return false;
  }

};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



struct Frustum {
  // silhouettes
  static const unsigned char edges[12][2];// = {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
  static const unsigned char faces[12][2];// = {{0,4},{0,2},{0,5},{0,3},{4,1},{2,1},{5,1},{3,1},{4,3},{2,4},{5,2},{3,5}};
  Vec3f vertices_[8];
  Vec3f planes_[6];
  float maxes_[6];
  Vec3f sil_[3][9]; // [x,y,z view dir][max 8 sil edges]
  Vec3f center_, lo_, hi_;
  float radius_;

  public:

  const Vec3f & vertex(const int i) const { return vertices_[i]; }
  const Vec3f & plane(const int i) const { return planes_[i]; }
  const Vec3f & center() const { return center_; }
  const float & radius() const { return radius_; }
  void makePlane(int, int, int, int);
  bool hitsAAB(const OBB &) const;

  Frustum();

  template<int c0, int c1>
  bool hitsAAB2D(const OBB & box) const {
    auto sils = sil_[(c0 + 2)%3];
    int pos = 0;
    while( sils[pos].x() != 0.0f || sils[pos].y() != 0.0f ) {
      const Vec3f & silv = sils[pos];
      float bmin = box.AA2DminimizeInDirection<c0,c1>(Vec2f(silv.x(), silv.y()));
      ++planeStatPerPair;
      if( silv.z() < bmin ) return false;
      ++pos;
    }
    return true;
  }
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


struct FullConvex {
  std::vector<Vec3f> vertices_;
  std::vector<Vec4us> edges_; // each stores: edgeV0, edgeV1, face0V, face1V
  struct Facet {
    Vec3f normal_;
    float maxAlong_;
    Facet(const Vec3f & n, float m) : normal_(n), maxAlong_(m) {}
  };
  std::vector<Facet> facets_;
  Vec3f center_;

  typedef std::map<const Polyhedron_3::Vertex_const_handle, int> VH_to_id;

  public:

  const Vec3f & vertex(const int i) const { return vertices_[i]; }
  const Vec3f & center() const { return center_; }
  size_t size() const { return vertices_.size(); }
  int nbVertices() const { return vertices_.size(); }
  int nbEdges() const { return edges_.size(); }
  int nbFacets() const { return facets_.size(); }

  FullConvex(int n, float shift) { exit(-1); }
  FullConvex(const Polyhedron_3 & ch);

  bool differenceCoversZeroInDir(const FullConvex & B, int & vA, int & vB, const Vec3f & dir) const {
    return false;
  }

  void differenceInDir(const FullConvex & B, int & vA, int & vB, float & maxOverA, float & minOverB, const Vec3f & dir) const {
    return;
  }

  bool differenceCoversZeroInDirLazy(const FullConvex & B, int & vA, int & vB, const Vec3f & dir) const {
    return false;
  }
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

struct VerticesOnly {
  std::vector<Vec3f> vertices_;
  Vec3f center_;

  public:

  const Vec3f & vertex(const int i) const { return vertices_[i]; }
  const Vec3f & center() const { return center_; }
  size_t size() const { return vertices_.size(); }
  size_t nbVertices() const { return vertices_.size(); }

  VerticesOnly(const Polyhedron_3 &) { exit(-1); }
  VerticesOnly(int n, float shift) {
    Vec3f s(shift, 0.0f, 0.0f);
    RandomPointOnSphere rps;
    Vec3f center;
    for( int i = 0; i < n; ++i ) {
      vertices_.push_back(rps.vec3f() + s);
      center = center + vertices_.back();
    }
    center_ = (1.0f / n) * center;
  }

  void makeFrustum() {
    vertices_.clear();
    vertices_.reserve(8);
    Frustum f;
    for( int i = 0; i < 8; ++i )
      vertices_.push_back(f.vertex(i));
    center_ = f.center();
  }

  void minimizeInDirection(const Vec3f & dir, float & mini, int & vA) const {
    mini = vertices_[0] | dir;
    vA = 0;
    const int na = vertices_.size();
    for( int i = 1; i < na; ++i ) {
      float tempA = vertices_[i] | dir;
      if( tempA < mini ) { mini = tempA; vA = i; }
    }
  }

  void maximizeInDirection(const Vec3f & dir, float & maxi, int & vA) const {
    maxi = vertices_[0] | dir;
    vA = 0;
    const int na = vertices_.size();
    for( int i = 1; i < na; ++i ) {
      float tempA = vertices_[i] | dir;
      if( tempA > maxi ) { maxi = tempA; vA = i; }
    }
  }

  void differenceInDir(const VerticesOnly & B, int & vA, int & vB, float & maxOverA, float & minOverB, const Vec3f & dir) const {
    // difference above is: A - B
    ++planeStatPerPair;
    maxOverA = vertices_[0] | dir;
    minOverB = B.vertices_[0] | dir;
    vA = vB = 0;
    const int na = vertices_.size();
    for( int i = 1; i < na; ++i ) {
      float tempA = vertices_[i] | dir;
      if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
    }
    const int nb = B.vertices_.size();
    for( int i = 1; i < nb; ++i ) {
      float tempB = B.vertices_[i] | dir;
      if( tempB < minOverB ) { minOverB = tempB; vB = i; }
    }
  }

  bool differenceCoversZeroInDir(const VerticesOnly & B, int & vA, int & vB, const Vec3f & dir) const {
    // difference above is: A - B
    ++planeStatPerPair;
    float maxOverA = vertices_[0] | dir;
    float minOverB = B.vertices_[0] | dir;
    vA = vB = 0;
    const int na = vertices_.size();
    for( int i = 1; i < na; ++i ) {
      float tempA = vertices_[i] | dir;
      if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
    }
    const int nb = B.vertices_.size();
    for( int i = 1; i < nb; ++i ) {
      float tempB = B.vertices_[i] | dir;
      if( tempB < minOverB ) { minOverB = tempB; vB = i; }
    }
    return maxOverA >= minOverB;
  }

  bool differenceCoversZeroInDirLazy(const VerticesOnly & B, int & vA, int & vB, const Vec3f & dir) const {
    // difference above is: A - B
    ++planeStatPerPair;
    float maxOverA = vertices_[0] | dir;
    float minOverB = B.vertices_[0] | dir;
    vA = vB = 0;
    //if( maxOverA >= minOverB ) return true;
    const int na = vertices_.size();
    const int nb = B.vertices_.size();
    int i;
    if( na <= nb ) {
      for( i = 1; i < na; ++i ) {
        float tempA = vertices_[i] | dir;
        float tempB = B.vertices_[i] | dir;
        if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
        if( tempB < minOverB ) { minOverB = tempB; vB = i; }
        if( maxOverA >= minOverB ) return true;
      }
      for( ; i < nb; ++i ) {
        float tempB = B.vertices_[i] | dir;
        if( tempB < minOverB ) {
          minOverB = tempB; vB = i;
          if( maxOverA >= minOverB ) return true;
        }
      }
    } else {
      for( i = 1; i < nb; ++i ) {
        float tempA = vertices_[i] | dir;
        float tempB = B.vertices_[i] | dir;
        if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
        if( tempB < minOverB ) { minOverB = tempB; vB = i; }
        if( maxOverA >= minOverB ) return true;
      }
      for( ; i < na; ++i ) {
        float tempA = vertices_[i] | dir;
        if( tempA > maxOverA ) {
          maxOverA = tempA; vA = i;
          if( maxOverA >= minOverB ) return true;
        }
      }
    }
    return false;
  }
};

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

struct SortedVertices {
  std::vector<Vec3f> vertices_;

  struct Node {
    Vec3f center_;
    float radius_;
    int split_;
    int firstChild_;
  };

  std::vector<Node> nodes_;

  Vec3f center_;

  public:

  const Vec3f & vertex(const int i) const { return vertices_[i]; }
  const Vec3f & center() const { return center_; }
  size_t size() const { return vertices_.size(); }
  size_t nbVertices() const { return vertices_.size(); }

  SortedVertices(const Polyhedron_3 &) { exit(-1); }
  SortedVertices(int nbv, float shift);
  void makeHierarchy();

  void maximizeInDirection(const Vec3f & dir, float & outMaxi, int & outV) const {
    float maxi = vertices_[0] | dir;
    int v = 0;
    struct TraverseData {
      int nodeId, start, end;
      TraverseData(int a, int b, int c) : nodeId(a), start(b), end(c) {}
    };
    const float dirLen = dir.length();
    static std::stack<TraverseData> stack;
    TraverseData td(0, 0, nbVertices());
    while( true ) {
      const Node & node = nodes_[td.nodeId];
      float maxOverSphere = (dir | node.center_) + node.radius_ * dirLen;
      if( maxOverSphere <= maxi ) goto get_work;
      if( node.firstChild_ == -1 ) { // no child
        for( int i = td.start; i < td.end; ++i ) {
          float temp = vertices_[i] | dir;
          if( temp > maxi ) {
            v = i;
            maxi = temp;
          }
        }
      } else {
        float t = dir | (nodes_[node.firstChild_].center_ - nodes_[node.firstChild_+1].center_);
        if( t > 0.0f ) {
          stack.emplace(node.firstChild_+1, node.split_, td.end);
          td = {node.firstChild_  , td.start, node.split_};
        } else {
          stack.emplace(node.firstChild_  , td.start, node.split_);
          td = {node.firstChild_+1, node.split_, td.end};
        }
        continue;
      }
get_work:
        if( stack.empty() ) { outV = v; outMaxi = maxi; return; }
        td = stack.top();
        stack.pop();
    }
  }

  bool differenceCoversZeroInDir(const SortedVertices & B, int & vA, int & vB, const Vec3f & dir) const {
    float maxOverA, minOverB;
    maximizeInDirection(dir, maxOverA, vA);
    Vec3f negdir(-dir.x(), -dir.y(), -dir.z());
    B.maximizeInDirection(negdir, minOverB, vB);
    ++planeStatPerPair;
    return maxOverA >= (-minOverB);
  }

  void differenceInDir(const SortedVertices & B, int & vA, int & vB, float & maxOverA, float & minOverB, const Vec3f & dir) const {
    maximizeInDirection(dir, maxOverA, vA);
    Vec3f negdir(-dir.x(), -dir.y(), -dir.z());
    B.maximizeInDirection(negdir, minOverB, vB);
    minOverB = - minOverB;
    ++planeStatPerPair;
  }

  bool differenceCoversZeroInDirLazy(const SortedVertices & B, int & vA, int & vB, const Vec3f & dir) const {
    return false;
  }
};

#endif // RAY_CONVEXES_H
