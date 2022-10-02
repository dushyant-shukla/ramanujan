#pragma once

#include "constants.h"

namespace voyager::math
{
/**
 * A 3D vector construct.
 */
struct vec3
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

    inline vec3() : x(0.0f), y(0.0f), z(0.0f) {}
    inline vec3(const float& _x, const float& _y, const float& _z) : x(_x), y(_y), z(_z) {}
    inline vec3(float* fv) : x(fv[0]), y(fv[1]), z(fv[3]) {}

    void normalize();
    bool operator==(const vec3& other);
    bool operator!=(const vec3& other);

    /**
     * Operators are defined as friend functions as we would like to keep the state of vectors unaffected
     * by the operations.
     */
    friend vec3  operator+(const vec3& a, const vec3& b);
    friend vec3  operator-(const vec3& a, const vec3& b);
    friend vec3  operator*(const vec3& a, float scaling_factor);
    friend vec3  operator*(float scaling_factor, const vec3& a);
    friend vec3  operator*(const vec3& a, const vec3& b);
    friend float dot(const vec3& a, const vec3& b);
    friend float lengthSq(const vec3& a);
    friend float length(const vec3& a);
    friend float distance(const vec3& a, const vec3& b);
    friend vec3  normalized(const vec3 a);

    /**
     * This method returns projection of vector A onto vector B.
     */
    friend vec3 projection(const vec3& a, const vec3& b);

    /**
     * This method returns rejection of vector A onto vector B.
     */
    friend vec3 rejection(const vec3& a, const vec3& b);

    /**
     * This function calculates the angle between two vectors in radians.
     */
    friend float angle(const vec3& a, const vec3& b);

    /**
     * This method returns reflection of vector A on a surface with
     * surface normal given by vector B.
     */
    friend vec3 reflection(const vec3& a, const vec3& b);

    friend vec3 cross(const vec3& a, const vec3& b);

    /**
     * .
     *
     * \param start
     * \param end
     * \param t
     * \return
     */
    friend vec3 lerp(const vec3& start, const vec3& end, float t);

    /**
     * .
     *
     * \param start
     * \param end
     * \param t
     * \return
     */
    friend vec3 slerp(const vec3& start, const vec3& end, float t);

    /**
     * .
     *
     * \param start
     * \param end
     * \param t
     * \return
     */
    friend vec3 nlerp(const vec3& start, const vec3& end, float t);
};
} // namespace voyager::math
