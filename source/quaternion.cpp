#include "quaternion.h"
#include "constants.h"

#include "quaternion.hpp"
#include <math.h>

namespace ramanujan
{

Quaternion& Quaternion::operator=(const Quaternion& other)
{
    this->scalar = other.scalar;
    this->vector = other.vector;
    return *this;

    // This implementation is causing troubles
    // Quaternion temp(other);
    // std::swap(*this, temp);
    // return *this;
}

Quaternion operator+(const Quaternion& a, const Quaternion& b)
{
    return Quaternion(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

Quaternion operator-(const Quaternion& a, const Quaternion& b)
{
    return Quaternion(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

Quaternion operator*(const Quaternion& q, real s)
{
    return Quaternion(q.x * s, q.y * s, q.z * s, q.w * s);
}

Quaternion operator-(const Quaternion& q)
{
    return Quaternion(-q.x, -q.y, -q.z, -q.w);
}

bool operator==(const Quaternion& left, const Quaternion& right)
{
    return (fabsf(left.x - right.x) <= Constants::EPSILON && fabsf(left.y - right.y) <= Constants::EPSILON &&
            fabsf(left.z - right.z) <= Constants::EPSILON && fabsf(left.w - right.w) <= Constants::EPSILON);
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

Quaternion operator^(const Quaternion& q, real t)
{
    real    angle   = 2.0f * acosf(q.scalar);
    Vector3 axis    = Normalized(q.vector);
    real    halfCos = cosf(t * angle * 0.5f);
    real    halfSin = sinf(t * angle * 0.5f);
    return Quaternion(axis.x * halfSin, axis.y * halfSin, axis.z * halfSin, halfCos);
}

bool SameOrientation(const Quaternion& left, const Quaternion& right)
{
    return (fabsf(left.x - right.x) <= Constants::EPSILON && fabsf(left.y - right.y) <= Constants::EPSILON &&
            fabsf(left.z - right.z) <= Constants::EPSILON && fabsf(left.w - right.w) <= Constants::EPSILON) ||
           (fabsf(left.x + right.x) <= Constants::EPSILON && fabsf(left.y + right.y) <= Constants::EPSILON &&
            fabsf(left.z + right.z) <= Constants::EPSILON && fabsf(left.w + right.w) <= Constants::EPSILON);
}

/**
 * The dot product measures how similar two quaternions are.
 */
real Dot(const Quaternion& a, const Quaternion& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

/**
 * The squared length of a quaternion is the same as the dot product of the quaternion with itself.
 */
real LengthSq(const Quaternion& q)
{
    return q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
}

real Length(const Quaternion& q)
{
    real length_sq = q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
    if(length_sq < Constants::EPSILON)
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
    real length_sq = q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
    if(length_sq < Constants::EPSILON)
    {
        return;
    }
    real inverse_length = 1.0f / sqrtf(length_sq);
    q.x *= inverse_length;
    q.y *= inverse_length;
    q.z *= inverse_length;
    q.w *= inverse_length;
}

Quaternion Normalized(const Quaternion& q)
{
    real length_sq = q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
    if(length_sq < Constants::EPSILON)
    {
        return Quaternion();
    }
    real inverse_length = 1.0f / sqrtf(length_sq);
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
    real length_sq = q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
    if(length_sq < Constants::EPSILON)
    {
        return Quaternion();
    }
    real reciprocal = 1.0f / length_sq;
    return Quaternion(-q.x * reciprocal, -q.y * reciprocal, -q.z * reciprocal, q.w * reciprocal);
}

Vector3 RotateVector(const Quaternion& q, const Vector3& v)
{
    return q * v;
}

Quaternion FromEulerAnglesRadians(real x, real y, real z)
{
    real c_x = cosf(x / 2.0f);
    real c_y = cosf(y / 2.0f);
    real c_z = cosf(z / 2.0f);

    real s_x = sinf(x / 2.0f);
    real s_y = sinf(y / 2.0f);
    real s_z = sinf(z / 2.0f);

    Quaternion q_x(s_x, 0, 0, c_x);
    Quaternion q_y(0, s_y, 0, c_y);
    Quaternion q_z(0, 0, s_z, c_z);

    return q_x * q_y * q_z; // Rotation order: XYZ
}

Quaternion FromEulerAnglesDegrees(real x, real y, real z)
{
    real rad_x = x * ramanujan::Constants::DEG_TO_RAD;
    real rad_y = y * ramanujan::Constants::DEG_TO_RAD;
    real rad_z = z * ramanujan::Constants::DEG_TO_RAD;
    return FromEulerAnglesRadians(rad_x, rad_y, rad_z);
}

Quaternion AngleAxis(real theta, const Vector3& axis)
{
    // Quaternions are created using an angle and an axis of rotation. A rotation of angle(theta) about an axis, can be
    // represented on a sphere as any directed arc whose length us 1/2 of angle(theta) on the plane perpendicular to the
    // rotation axis. Positive angles yield a counterclockwise rotation around the axis.
    Vector3 normalized_axis = Normalized(axis);

    // A quaternion can track two full rotations, which is 720 degrees. This makes the period of a quaternion 720
    // degrees. The period of sin/cos is 360 degrees. Dividing angle(theta) by 2 maps the range of a quaternion to the
    // range of sin/cos.
    real sin_theta = sinf(theta * 0.50f);

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

real GetAngle(const Quaternion& q)
{
    return 2.0f * acosf(q.w);
}

Matrix4 ToMatrix4(const Quaternion& q)
{
    Vector3 r = q * Vector3(1, 0, 0);
    Vector3 u = q * Vector3(0, 1, 0);
    Vector3 f = q * Vector3(0, 0, 1);
    return Matrix4(r.x, r.y, r.z, 0, u.x, u.y, u.z, 0, f.x, f.y, f.z, 0, 0, 0, 0, 1);
}

Quaternion ToQuaternion(const Matrix4& m)
{
    Vector3 up      = Normalized(Vector3(m.up.x, m.up.y, m.up.z));
    Vector3 forward = Normalized(Vector3(m.forward.x, m.forward.y, m.forward.z));
    Vector3 right   = Cross(up, forward);
    up              = Cross(forward, right);
    return LookRotation(forward, up);
}

Quaternion LookRotation(const Vector3& direction, const Vector3& up)
{
    // Find orthonormal basis vectors
    Vector3 f = Normalized(direction); // Object Forward
    Vector3 u = Normalized(up);        // Desired Up
    Vector3 r = Cross(u, f);           // Object Right
    u         = Cross(f, r);           // Object Up
    // From world forward to object forward
    Quaternion world_to_object = FromTo(Vector3(0, 0, 1), f);
    // what direction is the new object up?
    Vector3 object_up = world_to_object * Vector3(0, 1, 0);
    // From object up to desired up
    Quaternion u2u = FromTo(object_up, u);
    // Rotate to forward direction first
    // then twist to correct up
    Quaternion result = world_to_object * u2u;
    // Don't forget to normalize the result
    return Normalized(result);
}

Quaternion Mix(const Quaternion& from, const Quaternion& to, real t)
{
    return from * (1.0f - t) + to * t;
}

Quaternion Nlerp(const Quaternion& from, const Quaternion& to, real t)
{
    // Observe that Nlerp(from, to, t) = Normalized(Mix(from, to, t))
    return Normalized(from + (to - from) * t);
}

Quaternion Slerp(const Quaternion& from, const Quaternion& to, real t)
{
    if(fabsf(Dot(from, to)) > 1.0f - Constants::EPSILON)
    {
        return Nlerp(from, to, t);
    }
    Quaternion delta = Inverse(from) * to;
    return Normalized((delta ^ t) * from);
}

} // namespace ramanujan

namespace ramanujan::experimental
{

quat::quat() noexcept : x(real{0.0}), y(real{0.0}), z(real{0.0}), w(real{0.0}) {}

quat::quat(const real& _x, const real& _y, const real& _z, const real& _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}

quat::quat(const real* const data) noexcept : x(data[0]), y(data[1]), z(data[2]), w(data[3]) {}

quat::quat(const quat& other) noexcept : vector{other.vector}, scalar{other.scalar} {}

quat::quat(quat&& other) noexcept : vector{std::move(other.vector)}, scalar{std::move(other.scalar)} {}

quat& quat::operator=(const quat& other) noexcept
{
    if(*this != other)
    {
        x = other.x;
        y = other.y;
        z = other.z;
        w = other.w;
    }
    return *this;
}

//quat& quat::operator=(quat&& other) noexcept
//{
//    if(*this != other)
//    {
//        x = other.x;
//        y = other.y;
//        z = other.z;
//        w = other.w;
//    }
//    return *this;
//}

quat& quat::operator+=(const quat& rhs) noexcept
{
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    w += rhs.w;
    return *this;
}

quat& quat::operator-=(const quat& rhs) noexcept
{
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    w -= rhs.w;
    return *this;
}

quat& quat::operator*=(const quat& rhs) noexcept
{
    quat temp{*this};

    x = rhs.x * temp.w + rhs.y * temp.z - rhs.z * temp.y + rhs.w * temp.x;
    y = -rhs.x * temp.z + rhs.y * temp.w + rhs.z * temp.x + rhs.w * temp.y;
    z = rhs.x * temp.y - rhs.y * temp.x + rhs.z * temp.w + rhs.w * temp.z;
    w = -rhs.x * temp.x - rhs.y * temp.y - rhs.z * temp.z + rhs.w * temp.w;

    // scalar = temp.scalar * rhs.scalar - temp.vector.dot(rhs.vector);
    // vector = temp.scalar * rhs.vector + temp.vector * rhs.scalar + temp.vector.cross(rhs.vector);

    return *this;
}

quat& quat::operator/=(const quat& rhs) noexcept
{
    assert(rhs.x > kEpsilon);
    assert(rhs.y > kEpsilon);
    assert(rhs.z > kEpsilon);
    assert(rhs.w > kEpsilon);
    x /= rhs.x;
    y /= rhs.y;
    z /= rhs.z;
    w /= rhs.w;
    return *this;
}

quat quat::operator+(const quat& rhs) const noexcept
{
    return {x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w};
}

quat quat::operator-(const quat& rhs) const noexcept
{
    return {x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w};
}

quat quat::operator*(const quat& rhs) const noexcept
{
    quat result{};
    result.x = rhs.x * w + rhs.y * z - rhs.z * y + rhs.w * x;
    result.y = -rhs.x * z + rhs.y * w + rhs.z * x + rhs.w * y;
    result.z = rhs.x * y - rhs.y * x + rhs.z * w + rhs.w * z;
    result.w = -rhs.x * x - rhs.y * y - rhs.z * z + rhs.w * w;

    // result.scalar = scalar * rhs.scalar - vector.dot(rhs.vector);
    // result.vector = scalar * rhs.vector + vector * rhs.scalar + vector.cross(rhs.vector);

    return result;
}

quat quat::operator/(const quat& rhs) const noexcept
{
    assert(rhs.x > kEpsilon);
    assert(rhs.y > kEpsilon);
    assert(rhs.z > kEpsilon);
    assert(rhs.w > kEpsilon);
    return {x / rhs.x, y / rhs.y, z / rhs.z, w / rhs.w};
}

quat& quat::operator+=(const real& scalar) noexcept
{
    return *this;
}

quat& quat::operator-=(const real& scalar) noexcept
{
    return *this;
}

quat& quat::operator*=(const real& scalar) noexcept
{
    return *this;
}

quat& quat::operator/=(const real& scalar) noexcept
{
    return *this;
}

quat& quat::operator^=(const real& scalar) noexcept
{
    return *this;
}

quat quat::operator+(const real& scalar) const noexcept
{
    return quat();
}

quat quat::operator-(const real& scalar) const noexcept
{
    return quat();
}

quat quat::operator*(const real& scalar) const noexcept
{
    return quat();
}

quat quat::operator/(const real& scalar) const noexcept
{
    return quat();
}

quat quat::operator^(const real& scalar) const noexcept
{
    return quat();
}

quat quat::operator-() const noexcept
{
    return quat();
}

real* quat::data() noexcept
{
    return nullptr;
}

const real* quat::data() const noexcept
{
    return nullptr;
}

real quat::length() const noexcept
{
    return real();
}

real quat::lengthSquared() const noexcept
{
    return real();
}

quat& quat::normalize() noexcept
{
    return *this;
}

quat quat::normalized() const noexcept
{
    return quat();
}

bool quat::operator==(const quat& rhs) const noexcept
{
    return false;
}

bool quat::operator!=(const quat& rhs) const noexcept
{
    return false;
}

vec3 quat::operator*(const vec3& vector) const noexcept
{
    return vec3();
}

real quat::dot(const quat& rhs) const noexcept
{
    return real();
}

quat quat::conjugate() const noexcept
{
    return quat();
}

quat quat::inverse() const noexcept
{
    return quat();
}

vec3 quat::rotate(const vec3& vector) const noexcept
{
    return vec3();
}

vec3 quat::toEuler() const noexcept
{
    return vec3();
}

vec3 quat::toEulerDegrees() const noexcept
{
    return vec3();
}

quat quat::fromEuler(const vec3& euler) noexcept
{
    return quat();
}

quat quat::fromEulerDegrees(const vec3& euler) noexcept
{
    return quat();
}

quat quat::identity() noexcept
{
    return quat();
}

quat quat::lerp(const quat& a, const quat& b, const real& t) noexcept
{
    return quat();
}

quat quat::nlerp(const quat& a, const quat& b, const real& t) noexcept
{
    return quat();
}

quat quat::slerp(const quat& a, const quat& b, const real& t) noexcept
{
    return quat();
}

quat quat::fromAxisAngle(const vec3& axis, const real& angle) noexcept
{
    return quat();
}

quat quat::fromAxisAngle(const real& x, const real& y, const real& z, const real& angle) noexcept
{
    return quat();
}

quat quat::fromRotationMatrix(const mat4& matrix) noexcept
{
    return quat();
}

quat quat::fromDirection(const vec3& direction, const vec3& up) noexcept
{
    return quat();
}

quat quat::fromLookAt(const vec3& eye, const vec3& target, const vec3& up) noexcept
{
    return quat();
}

quat quat::fromRotationTo(const vec3& from, const vec3& to) noexcept
{
    return quat();
}

} // namespace ramanujan::experimental
