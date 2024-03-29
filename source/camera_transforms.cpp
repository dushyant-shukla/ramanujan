#include "camera_transforms.h"
#include "vector3.h"

#include <math.h>

namespace ramanujan
{

Matrix4 Frustrum(real left, real right, real bottom, real top, real near, real far)
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

Matrix4 Perspective(real fov, real aspect, real near, real far)
{
    real ymax = near * tanf(fov * 3.14159265359f / 360.0f);
    real xmax = ymax * aspect;
    return Frustrum(-xmax, xmax, -ymax, ymax, near, far);
}

Matrix4 Ortho(real left, real right, real bottom, real top, real near, real far)
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

Matrix4 LookAt(const Vector3& position, const Vector3& target, const Vector3& up)
{
    Vector3 f = Normalized(target - position) * -1.0f;
    Vector3 r = Cross(up, f); // Right handed
    if(r == Vector3(0, 0, 0))
    {
        return Matrix4(); // Error
    }
    Normalize(r);
    Vector3 u = Normalized(Cross(f, r)); // Right handed
    Vector3 t = Vector3(-Dot(r, position), -Dot(u, position), -Dot(f, position));
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