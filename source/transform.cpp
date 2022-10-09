#include "transform.h"
#include "constants.h"
#include "quaternion.h"
#include "vector3.h"

#include <math.h>

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
    Transform result;
    result.scale    = a.scale * b.scale;
    result.rotation = a.rotation * b.rotation;
    result.position = a.position * (a.scale * b.position);
    result.position = a.position + result.position;
    return result;
}

Transform Inverse(const Transform& t)
{
    Transform inverse;
    inverse.rotation = Inverse(t.rotation);
    inverse.scale.x  = fabsf(t.scale.x) < constants::EPSILON ? 0.0f : 1.0f / t.scale.x;
    inverse.scale.x  = fabsf(t.scale.y) < constants::EPSILON ? 0.0f : 1.0f / t.scale.y;
    inverse.scale.x  = fabsf(t.scale.z) < constants::EPSILON ? 0.0f : 1.0f / t.scale.z;

    Vector3 inverse_trans = t.position * (-1.0f);
    inverse.position      = inverse.rotation * (inverse.scale * inverse_trans);
    return inverse;
}

Transform Mix(const Transform& a, const Transform& b, float t)
{
    Quaternion b_rotation = b.rotation;
    if(Dot(a.rotation, b_rotation) < 0.0f)
    {
        b_rotation = -b_rotation;
    }
    return Transform(Lerp(a.position, b.position, t), Nlerp(a.rotation, b_rotation, t), Lerp(a.scale, b.scale, t));
}

// TODO
Matrix4 ToMatrix4(const Transform& t)
{
    // First, extract the rotation basis of the transform
    Vector3 x = t.rotation * Vector3(1, 0, 0);
    Vector3 y = t.rotation * Vector3(0, 1, 0);
    Vector3 z = t.rotation * Vector3(0, 0, 1);
    // Next, scale the basis vectors
    x = x * t.scale.x;
    y = y * t.scale.y;
    z = z * t.scale.z;
    // Extract the position of the transform
    Vector3 p = t.position;
    // Create matrix
    return Matrix4(x.x,
                   x.y,
                   x.z,
                   0, // X basis (& Scale)
                   y.x,
                   y.y,
                   y.z,
                   0, // Y basis (& scale)
                   z.x,
                   z.y,
                   z.z,
                   0, // Z basis (& scale)
                   p.x,
                   p.y,
                   p.z,
                   1 // Position
    );
}

// TODO
Transform ToTransform(const Matrix4& m)
{
    Transform out;
    out.position = Vector3(m.v[12], m.v[13], m.v[14]);
    out.rotation = ToQuaternion(m);
    Matrix4 rotScaleMat(m.v[0], m.v[1], m.v[2], 0, m.v[4], m.v[5], m.v[6], 0, m.v[8], m.v[9], m.v[10], 0, 0, 0, 0, 1);
    Matrix4 invRotMat    = ToMatrix4(Inverse(out.rotation));
    Matrix4 scaleSkewMat = rotScaleMat * invRotMat;
    out.scale            = Vector3(scaleSkewMat.v[0], scaleSkewMat.v[5], scaleSkewMat.v[10]);
    return out;
}

// TODO
Vector3 TransformPoint(const Transform& a, const Vector3& point)
{
    Vector3 out;
    out = a.rotation * (a.scale * point);
    out = a.position + out;
    return out;
}

// TODO
Vector3 TransformVector(const Transform& a, const Vector3& v)
{
    Vector3 out;
    out = a.rotation * (a.scale * v);
    return out;
}

} // namespace ramanujan