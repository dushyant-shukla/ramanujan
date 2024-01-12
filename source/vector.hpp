#ifndef RAMANUJAN_VECTOR
#define RAMANUJAN_VECTOR

#include "constants.h"
#include "precision.h"

#include <array>
#include <assert.h>
#include <ostream>

namespace ramanujan::experimental
{

template <typename VEC_TYPE, typename T, std::size_t N>
class TVector
{
public:
    using value_type      = T;
    using size_type       = std::size_t;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using pointer         = value_type*;
    using const_pointer   = const value_type*;

    VEC_TYPE&       type() { return static_cast<VEC_TYPE&>(*this); }
    const VEC_TYPE& type() const { return static_cast<const VEC_TYPE&>(*this); }

    VEC_TYPE& operator-() noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < N; ++i)
        {
            self.data[i] = -self.data[i];
        }
        return self;
    }

    VEC_TYPE& operator+=(const VEC_TYPE& rhs) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < N; ++i)
        {
            self.data[i] += rhs.data[i];
        }
        return self;
    }

    VEC_TYPE& operator-=(const VEC_TYPE& rhs) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < N; ++i)
        {
            self.data[i] -= rhs.data[i];
        }
        return self;
    }

    VEC_TYPE& operator*=(const VEC_TYPE& rhs) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < N; ++i)
        {
            self.data[i] *= rhs.data[i];
        }
        return self;
    }

    VEC_TYPE& operator/=(const VEC_TYPE& rhs) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < N; ++i)
        {
            self.data[i] /= rhs.data[i];
        }
        return self;
    }

    VEC_TYPE& operator+=(const T& scalar) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < N; ++i)
        {
            self.data[i] += scalar;
        }
        return self;
    }

    VEC_TYPE& operator-=(const T& scalar) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < N; ++i)
        {
            self.data[i] -= scalar;
        }
        return self;
    }

    VEC_TYPE& operator*=(const T& scalar) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < N; ++i)
        {
            self.data[i] *= scalar;
        }
        return self;
    }

    VEC_TYPE& operator/=(const T& scalar) noexcept
    {
        assert(scalar > T(0));
        auto& self = type();
        for(size_type i = 0; i < N; ++i)
        {
            self.data[i] /= scalar;
        }
        return self;
    }

    [[nodiscard]] VEC_TYPE operator+(const VEC_TYPE& rhs) const noexcept
    {
        auto&    self = type();
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = self.data[i] + rhs.data[i];
        }
        return result;
    }

    [[nodiscard]] VEC_TYPE operator-(const VEC_TYPE& rhs) const noexcept
    {
        auto&    self = type();
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = self.data[i] - rhs.data[i];
        }
        return result;
    }

    [[nodiscard]] VEC_TYPE operator*(const VEC_TYPE& rhs) const noexcept
    {
        auto&    self = type();
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = self.data[i] * rhs.data[i];
        }
        return result;
    }

    [[nodiscard]] VEC_TYPE operator/(const VEC_TYPE& rhs) const noexcept
    {
        auto&    self = type();
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = self.data[i] / rhs.data[i];
        }
        return result;
    }

    [[nodiscard]] VEC_TYPE operator+(const T& scalar) const noexcept
    {
        auto&    self = type();
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = self.data[i] + scalar;
        }
        return result;
    }

    [[nodiscard]] VEC_TYPE operator-(const T& scalar) const noexcept
    {
        auto&    self = type();
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = self.data[i] - scalar;
        }
        return result;
    }

    [[nodiscard]] VEC_TYPE operator*(const T& scalar) const noexcept
    {
        auto&    self = type();
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = self.data[i] * scalar;
        }
        return result;
    }

    [[nodiscard]] VEC_TYPE operator/(const T& scalar) const noexcept
    {
        assert(real_abs(scalar) > T(0));
        auto&    self = type();
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = self.data[i] / scalar;
        }
        return result;
    }

    [[nodiscard]] friend VEC_TYPE operator+(const T& scalar, const VEC_TYPE& rhs) noexcept
    {
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = rhs.data[i] + scalar;
        }
        return result;
    }

    [[nodiscard]] friend VEC_TYPE operator-(const T& scalar, const VEC_TYPE& rhs) noexcept
    {
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = rhs.data[i] - scalar;
        }
        return result;
    }

    [[nodiscard]] friend VEC_TYPE operator*(const T& scalar, const VEC_TYPE& rhs) noexcept
    {
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = rhs.data[i] * scalar;
        }
        return result;
    }

    [[nodiscard]] friend VEC_TYPE operator/(const T& scalar, const VEC_TYPE& rhs) noexcept
    {
        assert(scalar > T(0));
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = rhs.data[i] / scalar;
        }
        return result;
    }

    reference operator[](size_type i) noexcept
    {
        assert(i < N);
        return type().data[i];
    }

    const_reference operator[](size_type i) const noexcept
    {
        assert(i < N);
        return type().data[i];
    }

    void clear() noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < N; ++i)
        {
            self.data[i] = T(0);
        }
    }

    pointer data() noexcept { return type().data.data(); }

    const_pointer data() const noexcept { return type().data.data(); }

    [[nodiscard]] real length() const noexcept
    {
        real  result{0};
        auto& self = type();
        for(size_type i = 0; i < N; ++i)
        {
            result += self.data[i] * self.data[i];
        }
        return real_sqrt(result);
    }

    [[nodiscard]] real lengthSquared() const noexcept
    {
        real  result{0};
        auto& self = type();
        for(size_type i = 0; i < N; ++i)
        {
            result += self.data[i] * self.data[i];
        }
        return result;
    }

    [[nodiscard]] real magnitude() const noexcept { return length(); }

    [[nodiscard]] real magnitudeSquared() const noexcept { return lengthSquared(); }

    [[nodiscard]] bool operator==(const VEC_TYPE& rhs) const noexcept { return type().data == rhs.data; }

    [[nodiscard]] bool operator!=(const VEC_TYPE& rhs) const noexcept { return type().data != rhs.data; }

    [[nodiscard]] bool operator<(const VEC_TYPE& rhs) const noexcept { return type().data < rhs.data; }

    [[nodiscard]] bool operator<=(const VEC_TYPE& rhs) const noexcept { return type().data <= rhs.data; }

    [[nodiscard]] bool operator>(const VEC_TYPE& rhs) const noexcept { return type().data > rhs.data; }

    [[nodiscard]] bool operator>=(const VEC_TYPE& rhs) const noexcept { return type().data >= rhs.data; }

    VEC_TYPE& normalize() noexcept
    {
        auto& self = type();
        real  len  = length();
        if(len > 0.0f)
        {
            self /= len;
        }
        return self;
    }

    [[nodiscard]] VEC_TYPE normalized() const noexcept
    {
        VEC_TYPE result{};
        result.data = type().data;
        result.normalize();
        return result;
    }

    [[nodiscard]] value_type dot(const VEC_TYPE& rhs) const noexcept
    {
        auto&      self   = type();
        value_type result = value_type(0);
        for(size_type i = 0; i < N; ++i)
        {
            result += self.data[i] * rhs.data[i];
        }
        return result;
    }

    template <size_type D = N>
    [[nodiscard]] typename std::enable_if_t<D == 3, VEC_TYPE> cross(const VEC_TYPE& rhs) const noexcept
    {
        auto&    self = type();
        VEC_TYPE result{};
        result[0] = self.data[1] * rhs.data[2] - self.data[2] * rhs.data[1];
        result[1] = self.data[2] * rhs.data[0] - self.data[0] * rhs.data[2];
        result[2] = self.data[0] * rhs.data[1] - self.data[1] * rhs.data[0];
        return result;
    }

    [[nodiscard]] bool isZero() const noexcept { return type().lengthSquared() < kEpsilon; }

    [[nodiscard]] bool isParallel(const VEC_TYPE& b) const noexcept
    {
        return type().cross(b).lengthSquared() < kEpsilon;
    }

    [[nodiscard]] bool isOrthogonal(const VEC_TYPE& b) const noexcept { return type().dot(b) < kEpsilon; }

    /*!
     * @brief Linear interpolation can be calculated by scaling the difference between the two vectors,
     * and adding the result back to the original vector.
     *
     * The amount to lerp by is a normalized value (t) between 0 and 1. When t = 0, the interpolated
     * vector is the same as the starting vector. When t = 1, the interpolated vector is the same as
     * the end vector.
     *
     * Linearly interpolating between two vectors will always take the shortest path from one vector
     * to another.
     *
     * @param start The start vector
     * @param end The end vector
     * @param t The amount to lerp by
     * @return A linearly interpolated vector
     */
    [[nodiscard]] friend VEC_TYPE lerp(VEC_TYPE start, VEC_TYPE end, real t) noexcept
    {
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result[i] = start[i] + t * (end[i] - start[i]);
        }
        return result;
    }

    /*!
     * @brief Sometimes, the shortest path obtained between by linearly interpolating between two vectors
     * isn't the best path. Sometimes, we may want to interpolate between two vectors along the
     * shortest arc, i.e., Spherical Linear Interpolation (slerp).
     *
     * Assuming, theta is the angle between the two vectors, the formula for slerp is given by:
     *
     * slerp(start, end, t) = [sin((1-t)theta) / sin(theta) * start] + [sin((t)theta) / sin(theta) * end]
     *
     * @param start The start vector
     * @param end The end vector
     * @param t The amount to lerp by
     * @return A sphericaly interpolated vector
     */
    [[nodiscard]] friend VEC_TYPE slerp(VEC_TYPE start, VEC_TYPE end, real t) noexcept
    {
#ifdef COPILOT_GENERATED

        TVector result;
        real    dot = start.dot(end);
        if(dot < 0.0f)
        {
            end = -end;
            dot = -dot;
        }
        if(dot > 0.9995f)
        {
            return lerp(start, end, t);
        }
        real    theta           = acosf(dot) * t;
        TVector relative_vector = end - start * dot;
        relative_vector.normalize();
        return start * cosf(theta) + relative_vector * sinf(theta);

#endif

        /*
         * When the value of t is close to 0, slerp will yield unexpected results.
         * Therefore, we fall back to lerp or normalized lerp(nlerp).
         */
        if(t < real(0.01))
        {
            return lerp(start, end, t);
        }

        VEC_TYPE from = start.normalized();
        VEC_TYPE to   = end.normalized();

        real theta     = angle(from, to);
        real sin_theta = real_sin(theta);

        real a = real_sin((1.0f - t) * theta) / sin_theta;
        real b = real_sin(t * theta) / sin_theta;
        return a * from + b * to;
    }

    [[nodiscard]] friend VEC_TYPE nlerp(VEC_TYPE start, VEC_TYPE end, real t) noexcept
    {
        return lerp(start, end, t).normalized();
    }

    [[nodiscard]] friend real angle(const VEC_TYPE& a, const VEC_TYPE& b) noexcept
    {
        real sq_mag_a = a.lengthSquared();
        real sq_mag_b = b.lengthSquared();
        if(sq_mag_a < kEpsilon || sq_mag_b < kEpsilon)
        {
            return real(0.0);
        }
        return real_acos(a.dot(b) / (a.length() * b.length()));
    }

    [[nodiscard]] friend VEC_TYPE projection(const VEC_TYPE& a, const VEC_TYPE& b) noexcept
    {
        real sq_mag_b = b.lengthSquared();
        if(sq_mag_b < kEpsilon)
        {
            return VEC_TYPE(real(0.0));
        }

        /*
         * Basically, vector projection of vector 'a' on a non zero vector 'b' is given by [dot(a, unit(b)) * unit(b)].
         * Checkout https://en.wikipedia.org/wiki/Vector_projection
         */
        return b * (a.dot(b) / sq_mag_b);
    }

    [[nodiscard]] friend VEC_TYPE rejection(const VEC_TYPE& a, const VEC_TYPE& b) noexcept
    {
        /*
         * Rejection of vector 'a' onto vector 'b' is the opposite of projection of vector 'a' onto vector 'b'.
         * To find rejection of 'a' onto 'b', subtract the projection of 'a' onto 'b' from vector 'a'.
         */
        return a - projection(a, b);
    }

    [[nodiscard]] friend VEC_TYPE reflection(const VEC_TYPE& a, const VEC_TYPE& b) noexcept
    {
        return a - (2.0f * projection(a, b));
    }

    [[nodiscard]] friend bool orthonormalize(VEC_TYPE& a, VEC_TYPE& b, VEC_TYPE& c) noexcept
    {
#ifdef COPILOT_GENERATED
        a.normalize();
        b = rejection(b, a);
        b.normalize();
        c = rejection(rejection(c, a), b);
        c.normalize();
        return true;
#endif

        /*
         * Note that the construction of an orthonormal basis is a situation where it matters a great deal whether you
         * are working in a left - or right - handed coordinate system. The following algorithm is designed for right -
         * handed systems. If you need a left - handed coordinate system, then you can simply change the order of the
         * operands for both the cross - products. This will give you a left - handed orthonormal basis.
         */

#ifdef LEFT_HANDED_COORDINATE_SYSTEM

        // Left handed coordinate system
        a.normalize();
        c = b.cross(a);
        if(c.isZero())
        {
            return false;
        }
        c.normalize();
        b = a.cross(c);
        return true;

#else

        // Right handed coordinate system
        a.normalize();
        c = a.cross(b);
        if(c.isZero())
        {
            return false;
        }
        c.normalize();
        b = c.cross(a);
        return true;

#endif
    }

    friend std::ostream& operator<<(std::ostream& stream, const VEC_TYPE& vector)
    {
        stream << "[";
        for(size_type i = 0; i < N; ++i)
        {
            stream << vector.data[i];
            if(i < N - 1)
            {
                stream << ", ";
            }
        }
        stream << "]";
        return stream;
    }
};

struct vec2 : public TVector<vec2, real, 2>
{
    union
    {
        struct
        {
            real x, y;
        };
        std::array<real, 2> data;
    };

    vec2() noexcept : x(0), y(0)
    {
        data[0] = x;
        data[1] = y;
    }

    vec2(const real& v) noexcept : x(v), y(v)
    {
        data[0] = x;
        data[1] = y;
    }

    vec2(const real& _x, const real& _y) noexcept : x(_x), y(_y)
    {
        data[0] = x;
        data[1] = y;
    }

    vec2(const vec2& other) noexcept : x(other.x), y(other.y)
    {
        data[0] = x;
        data[1] = y;
    }

    vec2(vec2&& other) noexcept : x(other.x), y(other.y)
    {
        data[0] = x;
        data[1] = y;
    }

    vec2& operator=(const vec2& rhs) noexcept
    {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }

    vec2& operator=(vec2&& rhs) noexcept
    {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }
};

struct vec3 : public TVector<vec3, real, 3>
{
    union
    {
        struct
        {
            real x, y, z;
        };
        std::array<real, 3> data;
    };

    vec3() noexcept : x(0), y(0), z(0)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    vec3(const real& v) noexcept : x(v), y(v), z(v)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    vec3(const real& _x, const real& _y, const real& _z) noexcept : x(_x), y(_y), z(_z)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    vec3(const vec3& other) noexcept : x(other.x), y(other.y), z(other.z)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    vec3(vec3&& other) noexcept : x(other.x), y(other.y), z(other.z)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    vec3& operator=(const vec3& rhs) noexcept
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    vec3& operator=(vec3&& rhs) noexcept
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    // vec3& operator=(const vec3& rhs) noexcept = default;
    // vec3& operator=(vec3&& rhs) noexcept = default;

    [[nodiscard]] vec2 xy() const noexcept { return {x, y}; }

    [[nodiscard]] vec2 yz() const noexcept { return {y, z}; }

    [[nodiscard]] vec2 xz() const noexcept { return {x, z}; }
};

struct vec4 : public TVector<vec4, real, 4>
{
    union
    {
        struct
        {
            real x, y, z, w;
        };
        std::array<real, 4> data;
    };

    vec4() noexcept : x(0), y(0), z(0), w(0)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
        data[3] = w;
    }

    vec4(const real& v) noexcept : x(v), y(v), z(v), w(v)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
        data[3] = w;
    }

    vec4(const real& _x, const real& _y, const real& _z, const real& _w) noexcept : x(_x), y(_y), z(_z), w(_w)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
        data[3] = w;
    }

    vec4(const vec4& other) noexcept : x(other.x), y(other.y), z(other.z), w(other.w)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
        data[3] = w;
    }

    vec4(vec4&& other) noexcept : x(other.x), y(other.y), z(other.z), w(other.w)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
        data[3] = w;
    }

    vec4& operator=(const vec4& rhs) noexcept
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        w = rhs.w;
        return *this;
    }

    vec4& operator=(vec4&& rhs) noexcept
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        w = rhs.w;
        return *this;
    }

    vec4(const vec3& other, real w = real(1.0)) noexcept : x(other.x), y(other.y), z(other.z), w(w)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
        data[3] = w;
    }

    [[nodiscard]] vec2 xy() const noexcept { return {x, y}; }

    [[nodiscard]] vec2 yz() const noexcept { return {y, z}; }

    [[nodiscard]] vec2 xz() const noexcept { return {x, z}; }

    [[nodiscard]] vec3 xyz() const noexcept { return {x, y, z}; }

    [[nodiscard]] vec3 xzy() const noexcept { return {x, z, y}; }

    [[nodiscard]] vec3 yxz() const noexcept { return {y, x, z}; }

    [[nodiscard]] vec3 yzx() const noexcept { return {y, z, x}; }

    [[nodiscard]] vec3 zxy() const noexcept { return {z, x, y}; }

    [[nodiscard]] vec3 zyx() const noexcept { return {z, y, x}; }
};

struct ivec2 : public TVector<ivec2, int, 2>
{
    union
    {
        struct
        {
            int x, y;
        };
        std::array<int, 2> data;
    };

    ivec2() noexcept : x(0), y(0)
    {
        data[0] = x;
        data[1] = y;
    }

    ivec2(const unsigned& v) noexcept : x(v), y(v)
    {
        data[0] = x;
        data[1] = y;
    }

    ivec2(const unsigned& _x, const unsigned& _y) noexcept : x(_x), y(_y)
    {
        data[0] = x;
        data[1] = y;
    }

    ivec2(const ivec2& other) noexcept : x(other.x), y(other.y)
    {
        data[0] = x;
        data[1] = y;
    }

    ivec2(ivec2&& other) noexcept : x(other.x), y(other.y)
    {
        data[0] = x;
        data[1] = y;
    }

    ivec2& operator=(const ivec2& rhs) noexcept
    {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }

    ivec2& operator=(ivec2&& rhs) noexcept
    {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }
};

struct uvec2 : public TVector<uvec2, unsigned, 2>
{
    union
    {
        struct
        {
            unsigned x, y;
        };
        std::array<unsigned, 2> data;
    };

    uvec2() noexcept : x(0), y(0)
    {
        data[0] = x;
        data[1] = y;
    }

    uvec2(const unsigned& v) noexcept : x(v), y(v)
    {
        data[0] = x;
        data[1] = y;
    }

    uvec2(const unsigned& _x, const unsigned& _y) noexcept : x(_x), y(_y)
    {
        data[0] = x;
        data[1] = y;
    }

    uvec2(const uvec2& other) noexcept : x(other.x), y(other.y)
    {
        data[0] = x;
        data[1] = y;
    }

    uvec2(uvec2&& other) noexcept : x(other.x), y(other.y)
    {
        data[0] = x;
        data[1] = y;
    }

    uvec2& operator=(const uvec2& rhs) noexcept
    {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }

    uvec2& operator=(uvec2&& rhs) noexcept
    {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }
};

struct uvec3 : public TVector<uvec3, unsigned, 3>
{
    union
    {
        struct
        {
            unsigned x, y, z;
        };
        std::array<unsigned, 3> data;
    };

    uvec3() noexcept : x(0), y(0), z(0)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    uvec3(const unsigned& v) noexcept : x(v), y(v), z(v)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    uvec3(const unsigned& _x, const unsigned& _y, const unsigned& _z) noexcept : x(_x), y(_y), z(_z)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    uvec3(const uvec3& other) noexcept : x(other.x), y(other.y), z(other.z)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    uvec3(uvec3&& other) noexcept : x(other.x), y(other.y), z(other.z)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    uvec3& operator=(const uvec3& rhs) noexcept
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    uvec3& operator=(uvec3&& rhs) noexcept
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    // vec3& operator=(const vec3& rhs) noexcept = default;
    // vec3& operator=(vec3&& rhs) noexcept = default;

    [[nodiscard]] uvec2 xy() const noexcept { return {x, y}; }

    [[nodiscard]] uvec2 yz() const noexcept { return {y, z}; }

    [[nodiscard]] uvec2 xz() const noexcept { return {x, z}; }
};

} // namespace ramanujan::experimental

#endif