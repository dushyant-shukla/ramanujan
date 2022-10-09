#include "quaternion.h"
#include "constants.h"

#include <math.h>

namespace ramanujan
{

Quaternion operator+(const Quaternion& a, const Quaternion& b)
{
    return Quaternion(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

Quaternion operator-(const Quaternion& a, const Quaternion& b)
{
    return Quaternion(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

Quaternion operator*(const Quaternion& q, float s)
{
    return Quaternion(q.x * s, q.y * s, q.z * s, q.w * s);
}

Quaternion operator-(const Quaternion& q)
{
    return Quaternion(-q.x, -q.y, -q.z, -q.w);
}

bool operator==(const Quaternion& left, const Quaternion& right)
{
    return (fabsf(left.x - right.x) <= constants::EPSILON && fabsf(left.y - right.y) <= constants::EPSILON &&
            fabsf(left.z - right.z) <= constants::EPSILON && fabsf(left.w - right.w) <= constants::EPSILON);
}

bool operator!=(const Quaternion& left, const Quaternion& right)
{
    return !(left == right);
}

Quaternion operator*(const Quaternion& left, const Quaternion& right)
{
    return Quaternion(right.x * left.w + right.y * left.z - right.z * left.y + right.w * left.x,
                      -right.x * left.z + right.y * left.w + right.z * left.x + right.w * left.y,
                      right.x * left.y - right.y * left.x + right.z * left.w + right.w * left.z,
                      -right.x * left.x - right.y * left.y - right.z * left.z + right.w * left.w);
}

// Quaternion operator*(const Quaternion& left, const Quaternion& right)
//{
//     Quaternion result;
//     result.scalar = right.scalar * left.scalar - Dot(right.vector, left.vector);
//     result.vector = (left.vector * right.scalar) + (right.vector * left.scalar) + Cross(right.vector, left.vector);
//     return result;
// }

/**
 * Multiplying a vector by a quaternion will always yield a vector that is rotated by the quaternion.
 */
Vector3 operator*(const Quaternion& q, const Vector3& v)
{
    return q.vector * 2.0f * Dot(q.vector, v) + v * (q.scalar * q.scalar - Dot(q.vector, q.vector)) +
           Cross(q.vector, v) * 2.0f * q.scalar;
}

bool SameOrientation(const Quaternion& left, const Quaternion& right)
{
    return (fabsf(left.x - right.x) <= constants::EPSILON && fabsf(left.y - right.y) <= constants::EPSILON &&
            fabsf(left.z - right.z) <= constants::EPSILON && fabsf(left.w - right.w) <= constants::EPSILON) ||
           (fabsf(left.x + right.x) <= constants::EPSILON && fabsf(left.y + right.y) <= constants::EPSILON &&
            fabsf(left.z + right.z) <= constants::EPSILON && fabsf(left.w + right.w) <= constants::EPSILON);
}

/**
 * The dot product measures how similar two quaternions are.
 */
float Dot(const Quaternion& a, const Quaternion& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

/**
 * The squared length of a quaternion is the same as the dot product of the quaternion with itself.
 */
float LengthSq(const Quaternion& q)
{
    return q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
}

float Length(const Quaternion& q)
{
    float length_sq = q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
    if(length_sq < constants::EPSILON)
    {
        return 0.0f;
    }
    return sqrtf(length_sq);
}

/**
 * Quaternions that represent a rotatopm should always have a length of 1.
 */
void Normalize(Quaternion& q)
{
    float length_sq = q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
    if(length_sq < constants::EPSILON)
    {
        return;
    }
    float inverse_length = 1.0f / sqrtf(length_sq);
    q.x *= inverse_length;
    q.y *= inverse_length;
    q.z *= inverse_length;
    q.w *= inverse_length;
}

Quaternion Normalized(const Quaternion& q)
{
    float length_sq = q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
    if(length_sq < constants::EPSILON)
    {
        return Quaternion();
    }
    float inverse_length = 1.0f / sqrtf(length_sq);
    return Quaternion(q.x * inverse_length, q.y * inverse_length, q.z * inverse_length, q.w * inverse_length);
}

/**
 * The inverse of a normalized quaternion is its conjugate. The conjugate of a quaternion flips its axis of rotation.
 */
Quaternion Conjugate(const Quaternion& q)
{
    return Quaternion(-q.x, -q.y, -q.z, q.w);
}

/**
 * The inverse of a quaternion is the conjugate divided by the squared length of the quaternion.
 */
Quaternion Inverse(const Quaternion& q)
{
    float length_sq = q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
    if(length_sq < constants::EPSILON)
    {
        return Quaternion();
    }
    float reciprocal = 1.0f / length_sq;
    return Quaternion(-q.x * reciprocal, -q.y * reciprocal, -q.z * reciprocal, q.w * reciprocal);
}

Quaternion AngleAxis(float theta, const Vector3& axis)
{
    // Quaternions are created using an angle and an axis of rotation. A rotation of angle(theta) about an axis, can be
    // represented on a sphere as any directed arc whose length us 1/2 of angle(theta) on the plane perpendicular to the
    // rotation axis. Positive angles yield a counterclockwise rotation around the axis.
    Vector3 normalized_axis = Normalized(axis);

    // A quaternion can track two full rotations, which is 720 degrees. This makes the period of a quaternion 720
    // degrees. The period of sin/cos is 360 degrees. Dividing angle(theta) by 2 maps the range of a quaternion to the
    // range of sin/cos.
    float sin_theta = sinf(theta * 0.50f);

    return Quaternion(normalized_axis.x * sin_theta,
                      normalized_axis.y * sin_theta,
                      normalized_axis.z * sin_theta,
                      cosf(theta * 0.50f));
}

Quaternion FromTo(const Vector3& from, const Vector3& to)
{
    Vector3 f = Normalized(from);
    Vector3 t = Normalized(to);
    if(f == t)
    {
        return Quaternion();
    }
    else if(f == (t * (-1.0f)))
    {
        Vector3 ortho = Vector3(1.0f, 0, 0);
        if(fabsf(f.y) < fabsf(f.x))
        {
            ortho = Vector3(0, 1.0f, 0);
        }
        if(fabsf(f.z) < fabsf(f.y) && fabsf(f.z) < fabsf(f.x))
        {
            ortho = Vector3(0, 0, 1.0f);
        }
        Vector3 axis = Normalized(Cross(f, ortho));
        return Quaternion(axis.x, axis.y, axis.z, 0);
    }
    else
    {
        Vector3 half = Normalized(f + t);
        Vector3 axis = Cross(f, half);
        return Quaternion(axis.x, axis.y, axis.z, Dot(f, half));
    }
}

Vector3 GetAxis(const Quaternion& q)
{
    return Normalized(Vector3(q.x, q.y, q.z));
}

float GetAngle(const Quaternion& q)
{
    return 2.0f * acosf(q.w);
}

} // namespace ramanujan