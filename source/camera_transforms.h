#pragma once

#include "matrix4.h"

namespace ramanujan
{

struct Vector3;

Matrix4 Frustrum(float left, float right, float bottom, float top, float near, float far);
Matrix4 Perspective(float fov, float aspect, float near, float far);
Matrix4 Ortho(float left, float right, float bottom, float top, float near, float far);
Matrix4 LookAt(const Vector3& position, const Vector3& target, const Vector3& up);

} // namespace ramanujan