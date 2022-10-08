#include "transform.h"

namespace ramanujan
{

Transform::Transform(const Vector3& _position, const Quaternion& _rotation, const Vector3& _scale)
    : position(_position)
    , rotation(_rotation)
    , scale(_scale)
{
}

Transform Combine(const Transform& a, const Transform& b)
{
    return Transform();
}

Transform Inverse(const Transform& t)
{
    return Transform();
}

Transform Mix(const Transform& a, const Transform& b, float t)
{
    return Transform();
}

// TODO
Matrix4 ToMatrix4(const Transform& t)
{
    return Matrix4();
}

// TODO
Transform ToTransform(const Matrix4& m)
{
    return Transform();
}

// TODO
Vector3 TransformPoint(const Transform& a, const Vector3& point)
{
    return Vector3();
}

// TODO
Vector3 TransformVector(const Transform& a, const Vector3& v)
{
    return Vector3();
}

} // namespace ramanujan