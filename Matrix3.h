/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */

#ifndef RAY_MATRIX_3_H
#define RAY_MATRIX_3_H

#include "vec.h"
#include <functional>

template< typename NT >
struct GenMatrix3;

namespace Matrix {
  template< typename NT >
    void mul(const GenMatrix3<NT> & a, const GenMatrix3<NT> & b, GenMatrix3<NT> & r);
}

template< typename NT >
struct GenMatrix3 {
  NT d[3][3];

  typedef GenMatrix3<NT> Self;

  GenMatrix3(int) : d{
    {NT(1),NT(0),NT(0)},
    {NT(0),NT(1),NT(0)},
    {NT(0),NT(0),NT(1)}} {}

  GenMatrix3() : d{} {}

  GenMatrix3(const GenMatrix3 & m) {
    std::copy(m.d[0], m.d[0]+3, d[0]);
    std::copy(m.d[1], m.d[1]+3, d[1]);
    std::copy(m.d[2], m.d[2]+3, d[2]);
  }

  GenMatrix3(
      const NT i0, const NT i1, const NT i2,
      const NT i3, const NT i4, const NT i5,
      const NT i6, const NT i7, const NT i8)
  : d{{i0,i1,i2},{i3,i4,i5},{i6,i7,i8}} {}

  const NT & operator()(const int i,const int j) const {
    return d[i][j];
  }

  void col(int i, Vec3<NT> & v) const {
    v.set(d[0][i], d[1][i], d[2][i]);
  }

  NT det() const {
  return(
      d[0][0]*d[1][1]*d[2][2] + d[1][0]*d[2][1]*d[0][2] + d[2][0]*d[0][1]*d[1][2] -
      d[0][0]*d[2][1]*d[1][2] - d[1][0]*d[0][1]*d[2][2] - d[2][0]*d[1][1]*d[0][2]);
  }

  void inPlaceTranspose() {
    std::swap(d[0][1],d[1][0]);
    std::swap(d[0][2],d[2][0]);
    std::swap(d[1][2],d[2][1]);
  }

  void setRotation(const Vec3<NT> & inAxis, const NT angle) {
    NT cosa = ::cos(angle);
    NT sina = ::sin(angle);
    Vec3<NT> axis, u, v;
    axis = inAxis;
    axis.normalize();
    u = orthonormalVector(axis);
    v = Vec3<NT>::cross(axis, u);
    GenMatrix3<NT> localToGlobal(
        u.x(), v.x(), axis.x(),
        u.y(), v.y(), axis.y(),
        u.z(), v.z(), axis.z());
    GenMatrix3<NT> localRotation(
        cosa, -sina,  NT(0),
        sina,  cosa,  NT(0),
        NT(0), NT(0), NT(1));
    *this = localRotation;
    Matrix::mul(localToGlobal, (*this), localRotation);
    localToGlobal.inPlaceTranspose();
    Matrix::mul(localRotation, localToGlobal, (*this));
  }

  Vec3<NT> mulVec3(const Vec3<NT> & v) const {
    return Vec3<NT>(
        d[0][0] * v.x() + d[0][1] * v.y() + d[0][2] * v.z(),
        d[1][0] * v.x() + d[1][1] * v.y() + d[1][2] * v.z(),
        d[2][0] * v.x() + d[2][1] * v.y() + d[2][2] * v.z());
  }

  Vec3<NT> tMulVec3(const Vec3<NT> & v) const { // transposed multiplication
    return Vec3<NT>(
        d[0][0] * v.x() + d[1][0] * v.y() + d[2][0] * v.z(),
        d[0][1] * v.x() + d[1][1] * v.y() + d[2][1] * v.z(),
        d[0][2] * v.x() + d[1][2] * v.y() + d[2][2] * v.z());
  }

};

template< typename Out, typename NT >
Out & operator<<(Out & out, const GenMatrix3<NT> & m)
{
  for( int i = 0; i < 3; ++i ) {
    for( int j = 0; j < 3; ++j ) {
      out << m.d[i] << "\t\t";
    }
    out << '\n';
  }
  return out;
}


namespace Matrix {

  template< typename NT >
    void mul(const GenMatrix3<NT> & a, const GenMatrix3<NT> & b, GenMatrix3<NT> & __restrict r) {
      int pos(0);
      for( int i = 0; i < 3; ++i ) {
        for( int j = 0; j < 3; ++j )
          r.d[i][j] =
            a.d[i][0] * b.d[0][j] +
            a.d[i][1] * b.d[1][j] +
            a.d[i][2] * b.d[2][j];
      }
    }

  template< typename NT >
    void transMul(const GenMatrix3<NT> & a, const GenMatrix3<NT> & b, GenMatrix3<NT> & __restrict r) {
      int pos(0);
      for( int i = 0; i < 3; ++i ) {
        for( int j = 0; j < 3; ++j )
          r.d[i][j] =
            a.d[0][i] * b.d[0][j] +
            a.d[1][i] * b.d[1][j] +
            a.d[2][i] * b.d[2][j];
      }
    }

} // end of namespace Matrix

typedef GenMatrix3<float> Matrix3f;
typedef GenMatrix3<double> Matrix3d;

#endif // RAY_MATRIX_3_H
