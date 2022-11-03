#pragma once

#include "matrix4.h"
#include "vector3.h"

namespace ramanujan
{

/**
 * Quaternions are used to encode rotation about an arbitrary axis.
 */
struct Quaternion
{
    union
    {
        struct
        {
            float x;
            float y;
            float z;
            float w;
        };

        struct
        {
            Vector3 vector;
            float   scalar;
        };

        float v[4];
    };

    /**
     * Constructs an identity quaternion.
     */
    inline Quaternion() : x(0), y(0), z(0), w(1) {}

    inline Quaternion(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}

    inline Quaternion(float* fv) : x(fv[0]), y(fv[1]), z(fv[2]), w(fv[3]) {}
};

Quaternion operator+(const Quaternion& a, const Quaternion& b);
Quaternion operator-(const Quaternion& a, const Quaternion& b);
Quaternion operator*(const Quaternion& q, float s);
Quaternion operator-(const Quaternion& q);
bool       operator==(const Quaternion& left, const Quaternion& right);
bool       operator!=(const Quaternion& left, const Quaternion& right);
Quaternion operator*(const Quaternion& left, const Quaternion& right);
Vector3    operator*(const Quaternion& q, const Vector3& v);
Quaternion operator^(const Quaternion& q, float t);
bool       SameOrientation(const Quaternion& left, const Quaternion& right);
float      Dot(const Quaternion& a, const Quaternion& b);
float      LengthSq(const Quaternion& q);
float      Length(const Quaternion& q);
void       Normalize(Quaternion& q);
Quaternion Normalized(const Quaternion& q);
Quaternion Conjugate(const Quaternion& q);
Quaternion Inverse(const Quaternion& q);
Vector3    RotateVector(const Quaternion& q, const Vector3& v);

/**
 * This function constructs a quaternion from euler angles (in radians). The rotation order is XYZ.
 * @param x The angle on the x axis (in radians)
 * @param y The angle on the y axis (in radians)
 * @param z The angle on the z axis (in radians)
 */
Quaternion FromEulerAnglesRadians(float x, float y, float z);

/**
 * This function constructs a quaternion from euler angles (in degrees). The rotation order is XYZ.
 * @param x The angle on the x axis (in degrees)
 * @param y The angle on the y axis (in degrees)
 * @param z The angle on the z axis (in degrees)
 */
Quaternion FromEulerAnglesDegrees(float x, float y, float z);

/**
 * This method encodes an angle and an axis of rotation into a quaternion.
 *
 * \param theta The angle of rotation in radians
 * \param axis The axis of rotation
 * \return A quaternion representing the rotation
 */
Quaternion AngleAxis(float theta, const Vector3& axis);

/**
 * .
 */
Quaternion FromTo(const Vector3& from, const Vector3& to);

Vector3 GetAxis(const Quaternion& q);

float GetAngle(const Quaternion& q);

Matrix4 ToMatrix4(const Quaternion& q);

Quaternion ToQuaternion(const Matrix4& m);

Quaternion LookRotation(const Vector3& direction, const Vector3& up);

Quaternion Mix(const Quaternion& from, const Quaternion& to, float t);

Quaternion Nlerp(const Quaternion& from, const Quaternion& to, float t);

Quaternion Slerp(const Quaternion& from, const Quaternion& to, float t);

} // namespace ramanujan