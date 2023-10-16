#ifndef RAMANUJAN_VECTOR
#define RAMANUJAN_VECTOR

#include "constants.h"
#include "precision.h"

#include <array>

namespace ramanujan::experimental
{

template <typename T, std::size_t N>
class TVector
{
public:
    using value_type      = T;
    using size_type       = std::size_t;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using pointer         = value_type*;
    using const_pointer   = const value_type*;

    constexpr TVector() noexcept                          = default;
    constexpr TVector(const TVector&) noexcept            = default;
    constexpr TVector(TVector&&) noexcept                 = default;
    constexpr TVector& operator=(const TVector&) noexcept = default;
    constexpr TVector& operator=(TVector&&) noexcept      = default;

    template <std::size_t D = N>
    typename std::enable_if<D >= 4, TVector<T, 3>>::type xyz() const noexcept
    {
        return {m_data[0], m_data[1], m_data[2]};
    }

    template <std::size_t D = N>
    typename std::enable_if_t<D >= 3, TVector<T, 2>>::type xy() const noexcept
    {
        return {m_data[0], m_data[1]};
    }

    template <std::size_t D = N>
    typename std::enable_if_t<D >= 3, TVector<T, 2>>::type yz() const noexcept
    {
        return {m_data[1], m_data[2]};
    }

    template <std::size_t D = N>
    typename std::enable_if_t<D >= 3, TVector<T, 2>>::type xz() const noexcept
    {
        return {m_data[0], m_data[2]};
    }

    constexpr reference operator[](size_type index) noexcept { return m_data[index]; }

    constexpr const_reference operator[](size_type index) const noexcept { return m_data[index]; }

    constexpr reference at(size_type index) noexcept { return m_data.at(index); }

    constexpr const_reference at(size_type index) const noexcept { return m_data.at(index); }

    constexpr pointer data() noexcept { return m_data.data(); }

    constexpr const_pointer data() const noexcept { return m_data.data(); }

    constexpr size_type size() const noexcept { return m_data.size(); }

    constexpr bool empty() const noexcept { return m_data.empty(); }

    constexpr void fill(const T& value) noexcept { m_data.fill(value); }

    constexpr void swap(TVector& other) noexcept { m_data.swap(other.m_data); }

    constexpr TVector& operator+=(const TVector& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] += rhs.m_data[i];
        }
        return *this;
    }

    constexpr TVector& operator-=(const TVector& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] -= rhs.m_data[i];
        }
        return *this;
    }

    constexpr TVector& operator*=(const TVector& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] *= rhs.m_data[i];
        }
        return *this;
    }

    constexpr TVector& operator/=(const TVector& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] /= rhs.m_data[i];
        }
        return *this;
    }

    constexpr TVector& operator+=(const T& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] += rhs;
        }
        return *this;
    }

    constexpr TVector& operator-=(const T& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] -= rhs;
        }
        return *this;
    }

    constexpr TVector& operator*=(const T& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] *= rhs;
        }
        return *this;
    }

    constexpr TVector& operator/=(const T& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] /= rhs;
        }
        return *this;
    }

    constexpr TVector operator+(const TVector& rhs) const noexcept
    {
        TVector result(*this);
        result += rhs;
        return result;
    }

    constexpr TVector operator-(const TVector& rhs) const noexcept
    {
        TVector result(*this);
        result -= rhs;
        return result;
    }

    constexpr TVector operator*(const TVector& rhs) const noexcept
    {
        TVector result(*this);
        result *= rhs;
        return result;
    }

    constexpr TVector operator/(const TVector& rhs) const noexcept
    {
        TVector result(*this);
        result /= rhs;
        return result;
    }

    constexpr TVector operator+(const T& rhs) const noexcept
    {
        TVector result(*this);
        result += rhs;
        return result;
    }

    constexpr TVector operator-(const T& rhs) const noexcept
    {
        TVector result(*this);
        result -= rhs;
        return result;
    }

    constexpr TVector operator*(const T& rhs) const noexcept
    {
        TVector result(*this);
        result *= rhs;
        return result;
    }

    constexpr TVector operator/(const T& rhs) const noexcept
    {
        TVector result(*this);
        result /= rhs;
        return result;
    }

    constexpr bool operator==(const TVector& rhs) const noexcept { return m_data == rhs.m_data; }

    constexpr bool operator!=(const TVector& rhs) const noexcept { return m_data != rhs.m_data; }

    constexpr bool operator<(const TVector& rhs) const noexcept { return m_data < rhs.m_data; }

    constexpr bool operator<=(const TVector& rhs) const noexcept { return m_data <= rhs.m_data; }

    constexpr bool operator>(const TVector& rhs) const noexcept { return m_data > rhs.m_data; }

    constexpr bool operator>=(const TVector& rhs) const noexcept { return m_data >= rhs.m_data; }

    constexpr TVector& normalize() noexcept
    {
        T length = length();
        if(length > 0.0f)
        {
            *this /= length;
        }
        return *this;
    }

    constexpr TVector normalized() const noexcept
    {
        TVector result(*this);
        result.normalize();
        return result;
    }

    constexpr real length() const noexcept { return real_sqrt(lengthSquared()); }

    constexpr real lengthSquared() const noexcept { return *this.dot(*this); }

    constexpr real distance(const TVector& rhs) const noexcept { return (*this - rhs).length(); }

    constexpr real distanceSquared(const TVector& rhs) const noexcept { return (*this - rhs).lengthSquared(); }

    constexpr T dot(const TVector& rhs) const noexcept
    {
        T result = 0.0f;
        for(size_type i = 0; i < N; ++i)
        {
            result += m_data[i] * rhs.m_data[i];
        }
        return result;
    }

    template <std::size_t D = N>
    constexpr typename std::enable_if_t<D == 3, TVector<T, 3>> cross(const TVector& rhs) const noexcept
    {
        TVector result;
        result[0] = m_data[1] * rhs.m_data[2] - m_data[2] * rhs.m_data[1];
        result[1] = m_data[2] * rhs.m_data[0] - m_data[0] * rhs.m_data[2];
        result[2] = m_data[0] * rhs.m_data[1] - m_data[1] * rhs.m_data[0];
        return result;
    }

    constexpr bool isZero() const noexcept { return *this.lengthSquared() < kEpsilon; }

    constexpr bool isParallel(const TVector& b) const noexcept { return *this.cross(b).lengthSquared() < kEpsilon; }

    constexpr bool isOrthogonal(const TVector& b) const noexcept { return *this.dot(b) < kEpsilon; }

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
    friend constexpr TVector lerp(TVector start, TVector end, real t) noexcept
    {
        TVector result;
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
    friend constexpr TVector slerp(TVector start, TVector end, real t) noexcept
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
        if(t < 0.01f)
        {
            return lerp(start, end, t);
        }

        TVector from = start.normalized();
        TVector to   = end.normalized();

        real theta     = angle(from, to);
        real sin_theta = real_sin(theta);

        real a = real_sin((1.0f - t) * theta) / sin_theta;
        real b = real_sin(t * theta) / sin_theta;
        return a * from + b * to;
    }

    friend constexpr TVector nlerp(TVector start, TVector end, real t) noexcept
    {
        return lerp(start, end, t).normalized();
    }

    friend constexpr real angle(const TVector& a, const TVector& b) noexcept
    {
        real sq_mag_a = a.lengthSquared();
        real sq_mag_b = b.lengthSquared();
        if(sq_mag_a < kEpsilon || sq_mag_b < kEpsilon)
        {
            return 0.0f;
        }
        return real_acos(a.dot(b) / (a.length() * b.length()));
    }

    friend constexpr TVector projection(const TVector& a, const TVector& b) noexcept
    {
        real sq_mag_b = b.lengthSquared();
        if(sq_mag_b < kEpsilon)
        {
            return TVector(0.0f);
        }

        /*
         * Basically, vector projection of vector 'a' on a non zero vector 'b' is given by [dot(a, unit(b)) * unit(b)].
         * Checkout https://en.wikipedia.org/wiki/Vector_projection
         */
        return b * (a.dot(b) / sq_mag_b);
    }

    friend constexpr TVector rejection(const TVector& a, const TVector& b) noexcept
    {
        /*
         * Rejection of vector 'a' onto vector 'b' is the opposite of projection of vector 'a' onto vector 'b'.
         * To find rejection of 'a' onto 'b', subtract the projection of 'a' onto 'b' from vector 'a'.
         */
        return a - projection(a, b);
    }

    friend constexpr TVector reflection(const TVector& a, const TVector& b) noexcept
    {
        return a - (2.0f * projection(a, b));
    }

    friend constexpr bool orthonormalize(TVector& a, TVector& b, TVector& c) noexcept
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

#endif

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
    }

private:
    std::array<T, N> m_data;
};

template <typename T>
struct TVec2 : public TVector<T, 2>
{
    T x;
    T y;

    TVec2(const T& v) noexcept : x(v), y(v) {}
    TVec2(const T& _x, const T& _y) noexcept : x(_x), y(_y) {}
};

template <typename T>
struct TVec3 : public TVector<T, 3>
{
    T x;
    T y;
    T z;

    TVec3(const T& v) noexcept : x(v), y(v), z(v) {}
    TVec3(const T& _x, const T& _y, const T& _z) noexcept : x(_x), y(_y), z(_z) {}
};

template <typename T>
struct TVec4 : public TVector<T, 4>
{
    T x;
    T y;
    T z;
    T w;

    TVec4(const T& v) noexcept : x(v), y(v), z(v), w(v) {}
    TVec4(const T& _x, const T& _y, const T& _z, const T& _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
};

struct color3 : public TVector<float, 3>
{
    float r;
    float g;
    float b;

    color3(const float& v) noexcept : r(v), g(v), b(v) {}
    color3(const float& _r, const float& _g, const float& _b) noexcept : r(_r), g(_g), b(_b) {}
    color3(const float& _r, const float& _g, const float& _b, const float& _a) noexcept : r(_r), g(_g), b(_b) {}
};

struct color4 : public TVector<float, 4>
{
    float r;
    float g;
    float b;
    float a;

    color4(const float& v) noexcept : r(v), g(v), b(v), a(1.0f) {}
    color4(const float& _r, const float& _g, const float& _b) noexcept : r(_r), g(_g), b(_b), a(1.0f) {}
    color4(const float& _r, const float& _g, const float& _b, const float& _a) noexcept : r(_r), g(_g), b(_b), a(_a) {}
};

using vec2 = TVec2<real>;
using vec3 = TVec3<real>;
using vec4 = TVec4<real>;

using ivec2 = TVec2<int>;
using ivec3 = TVec3<int>;
using ivec4 = TVec4<int>;

using uvec2 = TVec2<unsigned>;
using uvec3 = TVec3<unsigned>;
using uvec4 = TVec4<unsigned>;

using bvec2 = TVec2<bool>;
using bvec3 = TVec3<bool>;
using bvec4 = TVec4<bool>;

} // namespace ramanujan::experimental

#endif