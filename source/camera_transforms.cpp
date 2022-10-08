#include "camera_transforms.h"
#include "vector3.h"

#include <math.h>

namespace ramanujan
{

Matrix4& CameraTransforms::Frustrum(float left, float right, float bottom, float top, float near, float far)
{
    if(left == right || top == bottom || near == far) // check for invalid boundary parameters
    {
        return Matrix4(); // return an identity matrix
    }
    return Matrix4((2.0f * near) / (right - left),
                   0,
                   0,
                   0,
                   0,
                   (2.0f * near) / (top - bottom),
                   0,
                   0,
                   (right + left) / (right - left),
                   (top + bottom) / (top - bottom),
                   (-(far + near)) / (far - near),
                   -1,
                   0,
                   0,
                   (-2 * far * near) / (far - near),
                   0);
}

Matrix4& CameraTransforms::Perspective(float fov, float aspect, float near, float far)
{
    float ymax = near * tanf(fov * 3.14159265359f / 360.0f);
    float xmax = ymax * aspect;
    return Frustrum(-xmax, xmax, -ymax, ymax, near, far);
}

Matrix4& CameraTransforms::Ortho(float left, float right, float bottom, float top, float near, float far)
{
    if(left == right || top == bottom || near == far) // check for invalid boundary conditions
    {
        return Matrix4(); // return an identity matrix
    }
    return Matrix4(2.0f / (right - left),
                   0,
                   0,
                   0,
                   0,
                   2.0f / (top - bottom),
                   0,
                   0,
                   0,
                   0,
                   -2.0f / (far - near),
                   0,
                   -((right + left) / (right - left)),
                   -((top + bottom) / (top - bottom)),
                   -((far + near) / (far - near)),
                   1);
}

Matrix4& CameraTransforms::LookAt(const Vector3& position, const Vector3& target, const Vector3& up)
{
    Vector3 f = normalized(target - position) * -1.0f;
    Vector3 r = cross(up, f); // Right handed
    if(r == Vector3(0, 0, 0))
    {
        return Matrix4(); // Error
    }
    r.normalize();
    Vector3 u = normalized(cross(f, r)); // Right handed
    Vector3 t = Vector3(-dot(r, position), -dot(u, position), -dot(f, position));
    return Matrix4(
        // Transpose upper 3x3 matrix to invert it
        r.x,
        u.x,
        f.x,
        0,
        r.y,
        u.y,
        f.y,
        0,
        r.z,
        u.z,
        f.z,
        0,
        t.x,
        t.y,
        t.z,
        1);
}

} // namespace ramanujan