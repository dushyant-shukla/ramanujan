#pragma once

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

    /**
     * This method encodes an angle and an axis of rotation into a quaternion.
     *
     * \param theta The angle of rotation
     * \param axis The axis of rotation
     * \return A quaternion representing the rotation
     */
    Quaternion AngleAxis(float theta, const Vector3& axis);
};

} // namespace ramanujan