/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */

#ifndef RAY_OBB_H
#define RAY_OBB_H

#include "vec.h"
#include "Matrix3.h"
#include "overlap.h"

struct alignas(16) OBB {
  Vec3f center_;
  Vec3f radii_;
  Matrix3f rotation_;
  float padding_;

  OBB() : center_(), radii_(), rotation_() {}

  OBB(const Vec3f & c, const Vec3f & r, const Matrix3f & m)
    : center_(c), radii_(r), rotation_(m) {}

  OBB(const Vec3f & c, const Vec3f & r)
    : center_(c), radii_(r) {
      Vec3f axis(drand48()-0.5f, drand48()-0.5f, drand48()-0.5f);
      rotation_.setRotation(axis, drand48() * M_PI);
    }

  const Vec3f & center() const { return center_; }
  Vec3f vertex() const { return center_ - rotation_.mulVec3(radii_); }
  Vec3f aaLo() const { return center_ - radii_; }
  Vec3f aaHi() const { return center_ + radii_; }

  void centeredAAmaximizeInDirection(const Vec3f & dir, float & maxi, Vec3f & maximizer) const {
    // ASSUMES center_ = (0,0,0)
    maximizer.set(
        dir.x() > 0.0f ? radii_.x() : - radii_.x(),
        dir.y() > 0.0f ? radii_.y() : - radii_.y(),
        dir.z() > 0.0f ? radii_.z() : - radii_.z()
        );
    maxi = dir | maximizer;
  }

  void AAmaximizeInDirection(const Vec3f & dir, float & maxi, Vec3f & maximizer) const {
    maximizer.set(
        center_.x() + (dir.x() > 0.0f ? radii_.x() : - radii_.x()),
        center_.y() + (dir.y() > 0.0f ? radii_.y() : - radii_.y()),
        center_.z() + (dir.z() > 0.0f ? radii_.z() : - radii_.z())
        );
    maxi = dir | maximizer;
  }

  void AAminimizeInDirection(const Vec3f & dir, float & mini, Vec3f & minimizer) const {
    minimizer.set(
        center_.x() + (dir.x() < 0.0f ? radii_.x() : - radii_.x()),
        center_.y() + (dir.y() < 0.0f ? radii_.y() : - radii_.y()),
        center_.z() + (dir.z() < 0.0f ? radii_.z() : - radii_.z())
        );
    mini = dir | minimizer;
  }

  template<int c0, int c1>
  float AA2DminimizeInDirection(const Vec2f & dir) const {
    Vec2f minimizer;
    minimizer.set(
        center_[c0] + (dir[0] < 0.0f ? radii_[c0] : - radii_[c0]),
        center_[c1] + (dir[1] < 0.0f ? radii_[c1] : - radii_[c1])
        );
    return (dir | minimizer);
  }

  void maximizeInDirection(const Vec3f & dir, float & maxi, Vec3f & maximizer) const {
    Vec3f localDir(rotation_.tMulVec3(dir));
    Vec3f diag(
        localDir.x() > 0.0f ? radii_.x() : - radii_.x(),
        localDir.y() > 0.0f ? radii_.y() : - radii_.y(),
        localDir.z() > 0.0f ? radii_.z() : - radii_.z()
        );
    maximizer = center_ + rotation_.mulVec3(diag);
    maxi = dir | maximizer;
  }

  void minimizeInDirection(const Vec3f & dir, float & mini, Vec3f & minimizer) const {
    Vec3f localDir(rotation_.tMulVec3(dir));
    Vec3f diag(
        localDir.x() > 0.0f ? radii_.x() : - radii_.x(),
        localDir.y() > 0.0f ? radii_.y() : - radii_.y(),
        localDir.z() > 0.0f ? radii_.z() : - radii_.z()
        );
    minimizer = center_ - rotation_.mulVec3(diag);
    mini = dir | minimizer;
  }
};

#endif // RAY_OBB_H
