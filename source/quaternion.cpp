#include "quaternion.h"

#include <cmath>

namespace ramanujan
{

Quaternion Quaternion::AngleAxis(float theta, const Vector3& axis)
{
    // Quaternions are created using an angle and an axis of rotation. A rotation of angle(theta) about an axis, can be
    // represented on a sphere as any directed arc whose length us 1/2 of angle(theta) on the plane perpendicular to the
    // rotation axis. Positive angles yield a counterclockwise rotation around the axis.
    Vector3 normalized_axis = normalized(axis);

    // A quaternion can track two full rotations, which is 720 degrees. This makes the period of a quaternion 720
    // degrees. The period of sin/cos is 360 degrees. Dividing angle(theta) by 2 maps the range of a quaternion to the
    // range of sin/cos.
    float sin_theta = sinf(theta * 0.50f);

    return Quaternion(normalized_axis.x * sin_theta,
                      normalized_axis.y * sin_theta,
                      normalized_axis.z * sin_theta,
                      cosf(theta * 0.50f));
}

} // namespace ramanujan