#pragma once

#include "constants.h"
#include "precision.h"
#include "vector.hpp"

namespace ramanujan
{
/*!
 * @brief A 3D vector construct.
 */
struct Vector3
{
    union
    {
        struct
        {
            real x;
            real y;
            real z;
        };
        real v[3];
    };

    inline Vector3() : x(0.0f), y(0.0f), z(0.0f) {}
    inline Vector3(real v) : x(v), y(v), z(v) {}
    inline Vector3(real _x, real _y, real _z) : x(_x), y(_y), z(_z) {}
    inline Vector3(real* fv) : x(fv[0]), y(fv[1]), z(fv[2]) {}

    bool           operator==(const Vector3& other);
    bool           operator!=(const Vector3& other);
    friend Vector3 operator+(const Vector3& a, const Vector3& b);
    friend Vector3 operator-(const Vector3& a, const Vector3& b);
    friend Vector3 operator*(const Vector3& a, real scaling_factor);
    friend Vector3 operator*(real scaling_factor, const Vector3& a);
    friend Vector3 operator*(const Vector3& a, const Vector3& b);

    /*!
     * @brief This method adds a scaled vector to another vector
     * @param a The vector to which the scaled vector is to be added
     * @param b The vector to be scaled
     * @param scaling_factor The value with which to scale the second vector
     * @return The vector resulting from adding a scaled vector 'b' to vector 'a'
     */
    friend Vector3 addScaledVector(const Vector3& a, const Vector3& b, const real& scaling_factor);

    /*!
     * @brief This method calculates and returns a component-wise product of two vectors.
     * @param a The first vector
     * @param b The second vector
     * @return The result of component-wise product between vector 'a' and vector 'b'
     */
    friend Vector3 componentProduct(const Vector3& a, const Vector3& b);

    operator ramanujan::experimental::vec3() const { return ramanujan::experimental::vec3{x, y, z}; }
};

real    Dot(const Vector3& a, const Vector3& b);
real    LengthSq(const Vector3& a);
real    Length(const Vector3& a);
real    Distance(const Vector3& a, const Vector3& b);
void    Normalize(Vector3& a);
Vector3 Normalized(const Vector3 a);

/*!
 * @brief This method returns projection of vector 'a' onto vector 'b'.
 * @param a The first vector
 * @param b The second vector
 * @return A vector representing the projection of vector 'a' onto vector 'b'
 */
Vector3 Projection(const Vector3& a, const Vector3& b);

/**
 * This method returns rejection of vector A onto vector B.
 */
Vector3 Rejection(const Vector3& a, const Vector3& b);

/**
 * This function calculates the angle between two vectors in radians.
 */
real Angle(const Vector3& a, const Vector3& b);

/**
 * This method returns reflection of vector A on a surface with
 * surface normal given by vector B.
 */
Vector3 Reflection(const Vector3& a, const Vector3& b);

Vector3 Cross(const Vector3& a, const Vector3& b);

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
Vector3 Lerp(const Vector3& start, const Vector3& end, real t);

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
Vector3 Slerp(const Vector3& start, const Vector3& end, real t);

/*!
 * @brief Normalized linear interpolation closely approximates spherical linear interpolation,
 * and is cheaper to calculate. Because of these reasons nlerp is much faster, and
 * generally, the better choice than slerp. Unlike slerp, nlerp is not constant in velocity. Therefore, use slerp if
 * constant interpolation velocity is required.
 *
 * @param start The start vector
 * @param end The end vector
 * @param t The amount to lerp by
 * @return A normalized linearly interpolated vector
 */
Vector3 Nlerp(const Vector3& start, const Vector3& end, real t);

bool makeOrthonormalBasis(Vector3& a, Vector3& b, Vector3& c);

} // namespace ramanujan
