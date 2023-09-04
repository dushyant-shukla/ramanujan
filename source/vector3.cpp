#include "vector3.h"

#include <math.h>

namespace ramanujan
{

bool Vector3::operator==(const Vector3& other)
{
    Vector3 diff(*this - other);
    return LengthSq(diff) < Constants::EPSILON;
}

bool Vector3::operator!=(const Vector3& other)
{
    return !(*this == other);
}

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

Vector3 operator*(const Vector3& a, real scaling_factor)
{
    Vector3 result(a.x * scaling_factor, a.y * scaling_factor, a.z * scaling_factor);
    return result;
}

Vector3 operator*(real scaling_factor, const Vector3& a)
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

Vector3 addScaledVector(const Vector3& a, const Vector3& b, const real& scaling_factor)
{
    return a + (scaling_factor * b);
}

Vector3 componentProduct(const Vector3& a, const Vector3& b)
{
    return a * b;
}

real Dot(const Vector3& a, const Vector3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

real LengthSq(const Vector3& a)
{
    return a.x * a.x + a.y * a.y + a.z * a.z;
}

real Length(const Vector3& a)
{
    real length_sq = a.x * a.x + a.y * a.y + a.z * a.z;
    if(length_sq < Constants::EPSILON)
    {
        return 0.0f;
    }
    return sqrtf(length_sq);
}

real Distance(const Vector3& a, const Vector3& b)
{
    return Length(a - b);
}

void Normalize(Vector3& a)
{
    real length_sq = a.x * a.x + a.y * a.y + a.z * a.z;
    if(length_sq < Constants::EPSILON)
    {
        return;
    }
    real inverted_length = 1.0f / sqrtf(length_sq);
    a.x *= inverted_length;
    a.y *= inverted_length;
    a.z *= inverted_length;
}

Vector3 Normalized(const Vector3 a)
{
    real length_sq = a.x * a.x + a.y * a.y + a.z * a.z;
    if(length_sq < Constants::EPSILON)
    {
        return a;
    }
    real inverted_length = 1.0f / sqrtf(length_sq);
    return Vector3(a.x * inverted_length, a.y * inverted_length, a.z * inverted_length);
}

Vector3 Projection(const Vector3& a, const Vector3& b)
{
    /*
     * To calculate the projection of vector 'a' onto vector 'b', vector 'a' must be broken down
     * into parallel and perpendicular components with respect to vector 'b'.
     *  1. Projection: The parallel component, i.e., length of vector 'a' in the direction of vector 'b'.
     *  2. Rejection: The perpendicular component, i.e., parallel component subtracted from vector 'a'.
     */

    real mag_b_sq = LengthSq(b);
    if(mag_b_sq < Constants::EPSILON)
    {
        return Vector3();
    }

    /*
     * If vector 'b' is a unit vector, then length of vector 'a' in the direction of vector 'b' is
     * given by a dot product between A and B. However, if neither input vector is normalized,
     * the dot product needs to be divided by the length of vector 'b' (the vector being projected onto).
     */
    real scale = Dot(a, b) / mag_b_sq;

    /*
     * Now that the parallel component of A with respect to B is known, vector 'b' can be scaled
     * by this component. Again, if B wasn't of unit length, the result will need to be divided
     * by the length of vector 'b' (the division by LengthSq(b) takes care of this).
     *
     * Basically, vector projection of vector 'a' on a non zero vector 'b' is given by [dot(a, unit(b)) * unit(b)].
     * Checkout https://en.wikipedia.org/wiki/Vector_projection
     */
    return b * scale;
}

Vector3 Rejection(const Vector3& a, const Vector3& b)
{
    /*
     * Rejection of vector 'a' onto vector 'b' is the opposite of projection of vector 'a' onto vector 'b'.
     * To find rejection of A onto B, subtract the projection of A onto B from vector 'a'.
     */
    Vector3 proj_a_on_b = Projection(a, b);
    return a - proj_a_on_b;
}

real Angle(const Vector3& a, const Vector3& b)
{
    real sq_mag_v1 = a.x * a.x + a.y * a.y + a.z * a.z;
    real sq_mag_v2 = b.x * b.x + b.y * b.y + b.z * b.z;

    if(sq_mag_v1 < Constants::EPSILON || sq_mag_v2 < Constants::EPSILON)
    {
        return 0.0f;
    }

    real dot    = a.x * b.x + a.y * b.y + a.z * b.z;
    real length = sqrtf(sq_mag_v1) * sqrtf(sq_mag_v2);
    return acosf(dot / length);
}

Vector3 Reflection(const Vector3& a, const Vector3& b)
{
    Vector3 proj_a_on_b = Projection(a, b);
    return a - (2 * proj_a_on_b);
}

Vector3 Cross(const Vector3& a, const Vector3& b)
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
Vector3 Lerp(const Vector3& start, const Vector3& end, real t)
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
Vector3 Slerp(const Vector3& start, const Vector3& end, real t)
{
    /*
     * When the value of t is close to 0, slerp will yield unexpected results.
     * Therefore, we fall back to lerp or normalized lerp(nlerp).
     */
    if(t < 0.01f)
    {
        return Lerp(start, end, t);
    }

    Vector3 from = Normalized(start);
    Vector3 to   = Normalized(end);

    real theta     = Angle(from, to);
    real sin_theta = sinf(theta);

    real a = sinf((1.0f - t) * theta) / sin_theta;
    real b = sinf(t * theta) / sin_theta;

    return a * from + b * to;
}

Vector3 Nlerp(const Vector3& start, const Vector3& end, real t)
{
    return Normalized(Lerp(start, end, t));
}

bool makeOrthonormalBasis(Vector3& a, Vector3& b, Vector3& c)
{
    /*
     * Note that the construction of an orthonormal basis is a situation where it matters a great deal whether you are
     * working in a left - or right - handed coordinate system. The following algorithm is designed for right - handed
     * systems. If you need a left - handed coordinate system, then you can simply change the order of the operands for
     * both the cross - products. This will give you a left - handed orthonormal basis.
     */

    Normalize(a);
    c             = Cross(a, b); // change to Cross(b, a) for a left-handed coordinate systems
    real c_mag_sq = LengthSq(c);
    if(c_mag_sq < Constants::EPSILON)
    {
        return false;
    }
    Normalize(c);
    b = Cross(c, a); // change to Cross(a, c) for a left-handed coordinate systems
    return true;
}

} // namespace ramanujan
