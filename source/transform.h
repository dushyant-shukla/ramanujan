#pragma once

#include "matrix4.h"
#include "quaternion.h"
#include "vector3.h"

namespace ramanujan
{

struct Transform
{
    Vector3    position;
    Quaternion rotation;
    Vector3    scale;

    inline Transform() : scale(Vector3(1.0f)) {}
    Transform(const Vector3& _position, const Quaternion& _rotation, const Vector3& _scale);
};

Transform Combine(const Transform& a, const Transform& b);
Transform Inverse(const Transform& t);
Transform Mix(const Transform& a, const Transform& b, real t);
Matrix4   ToMatrix4(const Transform& t);
Transform ToTransform(const Matrix4& m);
Vector3   TransformPoint(const Transform& a, const Vector3& point);
Vector3   TransformVector(const Transform& a, const Vector3& v);

} // namespace ramanujan