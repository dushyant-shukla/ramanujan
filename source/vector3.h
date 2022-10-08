#pragma once

#include "constants.h"

namespace ramanujan
{
/**
 * A 3D vector construct.
 */
struct Vector3
{
    union
    {
        struct
        {
            float x;
            float z;
            float y;
        };
        float v[3];
    };

    inline Vector3() : x(0.0f), y(0.0f), z(0.0f) {}
    inline Vector3(const float& _x, const float& _y, const float& _z) : x(_x), y(_y), z(_z) {}
    inline Vector3(float* fv) : x(fv[0]), y(fv[1]), z(fv[3]) {}

    bool operator==(const Vector3& other);
    bool operator!=(const Vector3& other);

    // Operators are defined as friend functions as we would like to keep the state of vectors unaffected by the
    // operations.

    friend Vector3 operator+(const Vector3& a, const Vector3& b);
    friend Vector3 operator-(const Vector3& a, const Vector3& b);
    friend Vector3 operator*(const Vector3& a, float scaling_factor);
    friend Vector3 operator*(float scaling_factor, const Vector3& a);
    friend Vector3 operator*(const Vector3& a, const Vector3& b);
    friend float   Dot(const Vector3& a, const Vector3& b);
    friend float   LengthSq(const Vector3& a);
    friend float   Length(const Vector3& a);
    friend float   Distance(const Vector3& a, const Vector3& b);
    friend void    Normalize(Vector3& a);
    friend Vector3 Normalized(const Vector3 a);

    /**
     * This method returns projection of vector A onto vector B.
     */
    friend Vector3 Projection(const Vector3& a, const Vector3& b);

    /**
     * This method returns rejection of vector A onto vector B.
     */
    friend Vector3 Rejection(const Vector3& a, const Vector3& b);

    /**
     * This function calculates the angle between two vectors in radians.
     */
    friend float Angle(const Vector3& a, const Vector3& b);

    /**
     * This method returns reflection of vector A on a surface with
     * surface normal given by vector B.
     */
    friend Vector3 Reflection(const Vector3& a, const Vector3& b);

    friend Vector3 Cross(const Vector3& a, const Vector3& b);

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
     * \param start The start vector
     * \param end The end vector
     * \param t The amount to lerp by
     * \return A linearly interpolated vector
     */
    friend Vector3 Lerp(const Vector3& start, const Vector3& end, float t);

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
    friend Vector3 Slerp(const Vector3& start, const Vector3& end, float t);

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
    friend Vector3 Nlerp(const Vector3& start, const Vector3& end, float t);
};
} // namespace ramanujan
