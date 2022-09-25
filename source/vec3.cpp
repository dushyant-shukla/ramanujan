#include "vec3.h"

namespace voyager::math
{
    vec3 operator+(const vec3& a, const vec3& b)
    {
        vec3 result(a.x + b.x, a.y + b.y, a.z + b.z);
        return result;
    }

    vec3 operator-(const vec3& a, const vec3& b)
    {
        vec3 result(a.x - b.x, a.y - b.y, a.z - b.z);
        return result;
    }

    vec3 operator*(const vec3& a, const float& scaling_factor)
    {
        vec3 result(a.x * scaling_factor, a.y * scaling_factor, a.z * scaling_factor);
        return result;
    }

    vec3 operator*(const vec3& a, const vec3& b)
    {
        vec3 result(a.x * b.x, a.y * b.y, a.z * b.z);
    }

    float dot(const vec3& a, const vec3& b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }
}