#include "vector3.h"

#include <math.h>

namespace ramanujan
{
Vector3 operator+(const Vector3& a, const Vector3& b)
{
    Vector3 result(a.x + b.x, a.y + b.y, a.z + b.z);
    return result;
}

Vector3 operator-(const Vector3& a, const Vector3& b)
{
    Vector3 result(a.x - b.x, a.y - b.y, a.z - b.z);
    return result;
}

Vector3 operator*(const Vector3& a, float scaling_factor)
{
    Vector3 result(a.x * scaling_factor, a.y * scaling_factor, a.z * scaling_factor);
    return result;
}

Vector3 operator*(float scaling_factor, const Vector3& a)
{
    Vector3 result(a.x * scaling_factor, a.y * scaling_factor, a.z * scaling_factor);
    return result;
}

/**
 * .
 *
 * \param a
 * \param b
 * \return
 */
Vector3 operator*(const Vector3& a, const Vector3& b)
{
    return Vector3(a.x * b.x, a.y * b.y, a.z * b.z);
}

/**
 * .
 *
 * \param a
 * \param b
 * \return
 */
float dot(const Vector3& a, const Vector3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/**
 * .
 *
 * \param a
 * \return
 */
float lengthSq(const Vector3& a)
{
    return a.x * a.x + a.y * a.y + a.z * a.z;
}

/**
 * .
 *
 * \param a
 * \return
 */
float length(const Vector3& a)
{
    float length_sq = a.x * a.x + a.y * a.y + a.z * a.z;
    if(length_sq < constants::EPSILON)
    {
        return 0.0f;
    }
    return sqrtf(length_sq);
}

/**
 * .
 *
 * \param a
 * \param b
 * \return
 */
float distance(const Vector3& a, const Vector3& b)
{
    return length(a - b);
}

/**
 * .
 *
 * \param a
 * \return
 */
Vector3 normalized(const Vector3 a)
{
    float length_sq = a.x * a.x + a.y * a.y + a.z * a.z;
    if(length_sq < constants::EPSILON)
    {
        return a;
    }
    float inverted_length = 1.0f / sqrtf(length_sq);
    return Vector3(a.x * inverted_length, a.y * inverted_length, a.z * inverted_length);
}

/**
 * To calculate the projection of vector A onto vector B, vector A must be broken down
 * into parallel and perpendicular components with respect to vector B.
 *  1. Projection: The parallel component, i.e., length of vector A in the direction of vector B.
 *  2. Rejection: The perpendicular component, i.e., paralle component subtracted from vector A.
 *
 * \param v The vector being projected onto.
 * \return The projection of current vector (*this) onto vector v.
 */
Vector3 projection(const Vector3& a, const Vector3& b)
{
    float mag_b_sq = lengthSq(b);
    if(mag_b_sq < constants::EPSILON)
    {
        return Vector3();
    }

    /*
     * If vector B is a unit vector, then length of vector A in the direction of vector B is
     * given by a dot product between A and B. However, if neither input vector is normalized,
     * the dot product needs to be divided by the length of vector B (the vector being projected onto).
     */
    float scale = dot(a, b) / mag_b_sq;

    /*
     * Now that the parallel component of A with respect to B is known, vector B can be scaled
     * by this component. Again, if B wasn't of unit length, the result will need to be divided
     * by the length of vector B.
     */
    return b * scale;
}

/**
 * Rejection of vector A onto vector B is the opposite of projection of vector A onto vector B.
 */
Vector3 rejection(const Vector3& a, const Vector3& b)
{
    /*
     * To find rejection of A onto B, subtract the projection of A onto B from vector A.
     */
    Vector3 proj_a_on_b = projection(a, b);
    return a - proj_a_on_b;
}

/**
 * .
 *
 * \param a
 * \param b
 * \return
 */
float angle(const Vector3& a, const Vector3& b)
{
    float sq_mag_v1 = a.x * a.x + a.y * a.y + a.z * a.z;
    float sq_mag_v2 = b.x * b.x + b.y * b.y + b.z * b.z;

    if(sq_mag_v1 < constants::EPSILON || sq_mag_v2 < constants::EPSILON)
    {
        return 0.0f;
    }

    float dot    = a.x * b.x + a.y * b.y + a.z * b.z;
    float length = sqrtf(sq_mag_v1) * sqrtf(sq_mag_v2);
    return acosf(dot / length);
}

/**
 * .
 *
 * \param a
 * \param b
 * \return
 */
Vector3 reflection(const Vector3& a, const Vector3& b)
{
    Vector3 proj_a_on_b = projection(a, b);
    return a - (2 * proj_a_on_b);
}

/**
 * .
 *
 * \param a
 * \param b
 * \return
 */
Vector3 cross(const Vector3& a, const Vector3& b)
{
    return Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

/**
 * Linear interpolation can be calculated by scaling the difference between the two vectors,
 * and adding the result back to the original vector.
 *
 * The amount to lerp by is a normalized value (t) between 0 and 1. When t = 0, the interpolated
 * vector is the same as the starting vector. When t = 1, the interpolated vector is the same as
 * the end vector.
 *
 * Linearly interpolating between two vectors will always take the shortest path from one vector
 * to another.
 *
 * \param a The start vector
 * \param b The end vector
 * \param t The amount to lerp by
 * \return A linearly interpolated vector
 */
Vector3 lerp(const Vector3& start, const Vector3& end, float t)
{
    return Vector3(start.x + (end.x - start.x) * t, start.y + (end.y - start.y) * t, start.z + (end.z - start.z) * t);
}

/**
 * Sometimes, the shortest path obtained between by linearly interpolating between two vectors
 * isn't the best path. Sometimes, we may want to interpolate between two vectors along the
 * shortest arc, i.e., Spherical Linear Interpolation (slerp).
 *
 * Assuming, theta is the angle between the two vectors, the formula for slerp is given by:
 *
 * slerp(start, end, t) = [sin((1-t)theta) / sin(theta) * start] + [sin((t)theta) / sin(theta) * end]
 *
 * \param start
 * \param end
 * \param t
 * \return A normalized linearly interpolated vector
 */
Vector3 slerp(const Vector3& start, const Vector3& end, float t)
{
    /*
     * When the value of t is close to 0, slerp will yield unexpected results.
     * Therefore, we fall back to lerp or normalized lerp(nlerp).
     */
    if(t < 0.01f)
    {
        return lerp(start, end, t);
    }

    Vector3 from = normalized(start);
    Vector3 to   = normalized(end);

    float theta     = angle(from, to);
    float sin_theta = sinf(theta);

    float a = sinf((1.0f - t) * theta) / sin_theta;
    float b = sinf(t * theta) / sin_theta;

    return a * from + b * to;
}

/**
 * Normalized linear interpolation closely approximates spherical linear interpolation,
 * and is cheaper to calculate. Because of these reasons nlerp is much faster, and
 * generally, the better choice than slerp.
 *
 * Unlike slerp, nlerp is not constant in velocity. Therefore, use slerp if constant
 * interpolation velocity is required.
 *
 * \param start
 * \param end
 * \param t
 * \return A normalized linearly interpolated vector
 */
Vector3 nlerp(const Vector3& start, const Vector3& end, float t)
{
    return normalized(lerp(start, end, t));
}

/**
 * .
 *
 * \return
 */
void Vector3::normalize()
{
    float length_sq = this->x * this->x + this->y * this->y + this->z * this->z;
    if(length_sq < constants::EPSILON)
    {
        return;
    }
    float inverted_length = 1.0f / sqrtf(length_sq);
    this->x *= inverted_length;
    this->y *= inverted_length;
    this->z *= inverted_length;
}

bool Vector3::operator==(const Vector3& other)
{
    Vector3 diff(*this - other);
    return lengthSq(diff) < constants::EPSILON;
}

bool Vector3::operator!=(const Vector3& other)
{
    return !(*this == other);
}
} // namespace ramanujan
