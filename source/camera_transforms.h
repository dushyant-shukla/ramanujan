#pragma once

#include "matrix4.h"

namespace ramanujan
{

struct Vector3;

Matrix4 Frustrum(real left, real right, real bottom, real top, real near, real far);
Matrix4 Perspective(real fov, real aspect, real near, real far);
Matrix4 Ortho(real left, real right, real bottom, real top, real near, real far);
Matrix4 LookAt(const Vector3& position, const Vector3& target, const Vector3& up);

} // namespace ramanujan