#ifndef RAMANUJAN_QUATERNION
#define RAMANUJAN_QUATERNION

#include "precision.h"
#include "vector.hpp"

#include <array>
#include <ostream>

namespace ramanujan::experimental
{

class quat
{

public:
    quat() noexcept;
    quat(const real& _x, const real& _y, const real& _z, const real& _w) noexcept;
    quat(const real* const _data) noexcept;
    quat(const quat& other) noexcept;
    quat(quat&& other) noexcept;

    quat& operator=(const quat& other) noexcept;
    //quat& operator=(quat&& other) noexcept;

    quat& operator+=(const quat& rhs) noexcept;
    quat& operator-=(const quat& rhs) noexcept;
    quat& operator*=(const quat& rhs) noexcept;
    quat& operator/=(const quat& rhs) noexcept;

    [[nodiscard]] quat operator+(const quat& rhs) const noexcept;
    [[nodiscard]] quat operator-(const quat& rhs) const noexcept;
    [[nodiscard]] quat operator*(const quat& rhs) const noexcept;
    [[nodiscard]] quat operator/(const quat& rhs) const noexcept;

    quat&              operator+=(const real& scalar) noexcept;
    quat&              operator-=(const real& scalar) noexcept;
    quat&              operator*=(const real& scalar) noexcept;
    quat&              operator/=(const real& scalar) noexcept;
    quat&              operator^=(const real& scalar) noexcept;
    [[nodiscard]] quat operator+(const real& scalar) const noexcept;
    [[nodiscard]] quat operator-(const real& scalar) const noexcept;
    [[nodiscard]] quat operator*(const real& scalar) const noexcept;
    [[nodiscard]] quat operator/(const real& scalar) const noexcept;
    [[nodiscard]] quat operator^(const real& scalar) const noexcept;

    [[nodiscard]] quat operator-() const noexcept;

    [[nodiscard]] real*       data() noexcept;
    [[nodiscard]] const real* data() const noexcept;

    [[nodiscard]] real length() const noexcept;
    [[nodiscard]] real lengthSquared() const noexcept;

    quat&              normalize() noexcept;
    [[nodiscard]] quat normalized() const noexcept;

    [[nodiscard]] bool operator==(const quat& rhs) const noexcept;
    [[nodiscard]] bool operator!=(const quat& rhs) const noexcept;

    [[nodiscard]] vec3 operator*(const vec3& vector) const noexcept;

    [[nodiscard]] real dot(const quat& rhs) const noexcept;

    [[nodiscard]] quat conjugate() const noexcept;
    [[nodiscard]] quat inverse() const noexcept;

    [[nodiscard]] vec3 rotate(const vec3& vector) const noexcept;

    [[nodiscard]] vec3 toEuler() const noexcept;
    [[nodiscard]] vec3 toEulerDegrees() const noexcept;

    [[nodiscard]] static quat fromEuler(const vec3& euler) noexcept;
    [[nodiscard]] static quat fromEulerDegrees(const vec3& euler) noexcept;

    [[nodiscard]] static quat identity() noexcept;

    [[nodiscard]] static quat lerp(const quat& a, const quat& b, const real& t) noexcept;
    [[nodiscard]] static quat nlerp(const quat& a, const quat& b, const real& t) noexcept;
    [[nodiscard]] static quat slerp(const quat& a, const quat& b, const real& t) noexcept;

    [[nodiscard]] static quat fromAxisAngle(const vec3& axis, const real& angle) noexcept;
    [[nodiscard]] static quat fromAxisAngle(const real& x, const real& y, const real& z, const real& angle) noexcept;

    [[nodiscard]] static quat fromRotationMatrix(const mat4& matrix) noexcept;

    [[nodiscard]] static quat fromDirection(const vec3& direction, const vec3& up) noexcept;

    [[nodiscard]] static quat fromLookAt(const vec3& eye, const vec3& target, const vec3& up) noexcept;

    [[nodiscard]] static quat fromRotationTo(const vec3& from, const vec3& to) noexcept;

public:
    union
    {
        struct
        {
            real x, y, z, w;
        };

        struct
        {
            vec3 vector;
            real scalar;
        };
    };
};

} // namespace ramanujan::experimental

#endif //  RAMANUJAN_QUATERNION
