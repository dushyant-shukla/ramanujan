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
            self.data[i] *= scalar;
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
        assert(scalar > T(0));
        auto&    self = type();
        VEC_TYPE result{};
        for(size_type i = 0; i < N; ++i)
        {
            result.data[i] = self.data[i] / scalar;
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
        result.data = data;
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
    [[nodiscard]] typename std::enable_if_t<D == 3, VEC_TYPE> cross(const TVector& rhs) const noexcept
    {
        auto&    self = type();
        VEC_TYPE result{};
        result[0] = self.data[1] * rhs.m_data[2] - self.data[2] * rhs.m_data[1];
        result[1] = self.data[2] * rhs.m_data[0] - self.data[0] * rhs.m_data[2];
        result[2] = self.data[0] * rhs.m_data[1] - self.data[1] * rhs.m_data[0];
        return result;
    }

    [[nodiscard]] bool isZero() const noexcept { return type().lengthSquared() < kEpsilon; }

    [[nodiscard]] bool isParallel(const VEC_TYPE& b) const noexcept
    {
        return type().cross(b).lengthSquared() < kEpsilon;
    }

    [[nodiscard]] bool isOrthogonal(const VEC_TYPE& b) const noexcept { return type().dot(b) < kEpsilon; }

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

    //vec3& operator=(const vec3& rhs) noexcept = default;
    //vec3& operator=(vec3&& rhs) noexcept = default;

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

} // namespace ramanujan::experimental

#endif