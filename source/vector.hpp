#ifndef RAMANUJAN_VECTOR
#define RAMANUJAN_VECTOR

#include "constants.h"
#include "precision.h"

#include <array>

namespace ramanujan::experimental
{

template <typename TYPE, typename T, std::size_t N>
class TVector
{
public:
    using value_type      = T;
    using size_type       = std::size_t;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using pointer         = value_type*;
    using const_pointer   = const value_type*;

    TVector() noexcept               = default;
    TVector(const TVector&) noexcept = default;
    TVector(TVector&&) noexcept      = default;

    TVector(const T& value) noexcept { m_data.fill(value); }

    template <typename... TArgs>
    TVector(TArgs... args) noexcept
    {
        static_assert(N == sizeof...(args));
        int j = 0;
        for(auto value : std::initializer_list<std::common_type_t<TArgs...>>{args...})
        {
            m_data[j] = value;
            ++j;
        }
    }

    TVector& operator=(const TVector& rhs) noexcept
    {
        if(this != &rhs)
        {
            m_data = rhs.m_data;
        }
        return **static_cast<TYPE*>(this);
    }

    TVector& operator=(TVector&& rhs) noexcept
    {
        if(this != &rhs)
        {
            m_data = rhs.m_data;
        }
        return *static_cast<TYPE*>(this);
    }

    template <std::size_t D = N>
    typename std::enable_if<D >= 4, TVector<TYPE, T, 3>>::type xyz() const noexcept
    {
        return {m_data[0], m_data[1], m_data[2]};
    }

    template <std::size_t D = N>
    typename std::enable_if_t<D >= 3, TVector<TYPE, T, 2>>::type xy() const noexcept
    {
        return {m_data[0], m_data[1]};
    }

    template <std::size_t D = N>
    typename std::enable_if_t<D >= 3, TVector<TYPE, T, 2>>::type yz() const noexcept
    {
        return {m_data[1], m_data[2]};
    }

    template <std::size_t D = N>
    typename std::enable_if_t<D >= 3, TVector<TYPE, T, 2>>::type xz() const noexcept
    {
        return {m_data[0], m_data[2]};
    }

    reference operator[](size_type index) noexcept { return m_data[index]; }

    const_reference operator[](size_type index) const noexcept { return m_data[index]; }

    reference at(size_type index) noexcept { return m_data.at(index); }

    const_reference at(size_type index) const noexcept { return m_data.at(index); }

    pointer data() noexcept { return m_data.data(); }

    const_pointer data() const noexcept { return m_data.data(); }

    size_type size() const noexcept { return m_data.size(); }

    bool empty() const noexcept { return m_data.empty(); }

    void fill(const T& value) noexcept { m_data.fill(value); }

    void swap(TVector& other) noexcept { m_data.swap(other.m_data); }

    void clear() noexcept { m_data.fill(0.0f); }

    TVector& operator+=(const TVector& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] += rhs.m_data[i];
        }
        return *static_cast<TYPE*>(this);
    }

    TVector& operator-=(const TVector& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] -= rhs.m_data[i];
        }
        return *static_cast<TYPE*>(this);
    }

    TVector& operator*=(const TVector& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] *= rhs.m_data[i];
        }
        return *static_cast<TYPE*>(this);
    }

    TVector& operator/=(const TVector& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] /= rhs.m_data[i];
        }
        return *static_cast<TYPE*>(this);
    }

    TVector& operator+=(const T& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] += rhs;
        }
        return *static_cast<TYPE*>(this);
    }

    TVector& operator-=(const T& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] -= rhs;
        }
        return *static_cast<TYPE*>(this);
    }

    TVector& operator*=(const T& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] *= rhs;
        }
        return *static_cast<TYPE*>(this);
    }

    TVector& operator/=(const T& rhs) noexcept
    {
        for(size_type i = 0; i < N; ++i)
        {
            m_data[i] /= rhs;
        }
        return *static_cast<TYPE*>(this);
    }

    [[nodiscard]] friend TYPE operator+(const TYPE& lhs, const TYPE& rhs) noexcept
    {
        TYPE result = lhs;
        result += rhs;
        return result;
    }

    [[nodiscard]] friend TYPE operator-(const TYPE& lhs, const TYPE& rhs) noexcept
    {
        TYPE result = lhs;
        result -= rhs;
        return result;
    }

    [[nodiscard]] friend TYPE operator*(const TYPE& lhs, const TYPE& rhs) noexcept
    {
        TYPE result = lhs;
        result *= rhs;
        return result;
    }

    [[nodiscard]] friend TYPE operator/(const TYPE& lhs, const TYPE& rhs) noexcept
    {
        TYPE result = lhs;
        result /= rhs;
        return result;
    }

    [[nodiscard]] friend TYPE operator+(const TYPE& lhs, const T& rhs) noexcept
    {
        TYPE result = lhs;
        result += rhs;
        return result;
    }

    [[nodiscard]] friend TYPE operator-(const TYPE& lhs, const T& rhs) noexcept
    {
        TYPE result = lhs;
        result -= rhs;
        return result;
    }

    [[nodiscard]] friend TYPE operator*(const TYPE& lhs, const T& rhs) noexcept
    {
        TYPE result = lhs;
        result *= rhs;
        return result;
    }

    [[nodiscard]] friend TYPE operator/(const TYPE& lhs, const T& rhs) noexcept
    {
        TYPE result = lhs;
        result /= rhs;
        return result;
    }

    [[nodiscard]] bool operator==(const TVector& rhs) const noexcept { return m_data == rhs.m_data; }

    [[nodiscard]] bool operator!=(const TVector& rhs) const noexcept { return m_data != rhs.m_data; }

    [[nodiscard]] bool operator<(const TVector& rhs) const noexcept { return m_data < rhs.m_data; }

    [[nodiscard]] bool operator<=(const TVector& rhs) const noexcept { return m_data <= rhs.m_data; }

    [[nodiscard]] bool operator>(const TVector& rhs) const noexcept { return m_data > rhs.m_data; }

    [[nodiscard]] bool operator>=(const TVector& rhs) const noexcept { return m_data >= rhs.m_data; }

    TVector& normalize() noexcept
    {
        T len = length();
        if(len > 0.0f)
        {
            *this /= len;
        }
        return *static_cast<TYPE*>(this);
    }

    [[nodiscard]] TYPE normalized() const noexcept
    {
        TYPE result{};
        result.m_data = m_data;
        result.normalize();
        return result;
    }

    [[nodiscard]] real length() const noexcept { return real_sqrt(lengthSquared()); }

    [[nodiscard]] real lengthSquared() const noexcept { return dot(*this); }

    [[nodiscard]] real distance(const TVector& rhs) const noexcept { return (*this - rhs).length(); }

    [[nodiscard]] real distanceSquared(const TVector& rhs) const noexcept { return (*this - rhs).lengthSquared(); }

    [[nodiscard]] T dot(const TVector& rhs) const noexcept
    {
        T result = 0.0f;
        for(size_type i = 0; i < N; ++i)
        {
            result += m_data[i] * rhs.m_data[i];
        }
        return result;
    }

    template <std::size_t D = N>
    [[nodiscard]] typename std::enable_if_t<D == 3, TYPE> cross(const TVector& rhs) const noexcept
    {
        TYPE result;
        result[0] = m_data[1] * rhs.m_data[2] - m_data[2] * rhs.m_data[1];
        result[1] = m_data[2] * rhs.m_data[0] - m_data[0] * rhs.m_data[2];
        result[2] = m_data[0] * rhs.m_data[1] - m_data[1] * rhs.m_data[0];
        return result;
    }

    [[nodiscard]] bool isZero() const noexcept { return lengthSquared() < kEpsilon; }

    [[nodiscard]] bool isParallel(const TYPE& b) const noexcept { return cross(b).lengthSquared() < kEpsilon; }

    [[nodiscard]] bool isOrthogonal(const TYPE& b) const noexcept { return dot(b) < kEpsilon; }

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
    [[nodiscard]] friend TYPE lerp(TYPE start, TYPE end, real t) noexcept
    {
        TYPE result;
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
    [[nodiscard]] friend TYPE slerp(TYPE start, TYPE end, real t) noexcept
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

        TYPE from = start.normalized();
        TYPE to   = end.normalized();

        real theta     = angle(from, to);
        real sin_theta = real_sin(theta);

        real a = real_sin((1.0f - t) * theta) / sin_theta;
        real b = real_sin(t * theta) / sin_theta;
        return a * from + b * to;
    }

    [[nodiscard]] friend TYPE nlerp(TYPE start, TYPE end, real t) noexcept { return lerp(start, end, t).normalized(); }

    [[nodiscard]] friend real angle(const TYPE& a, const TYPE& b) noexcept
    {
        real sq_mag_a = a.lengthSquared();
        real sq_mag_b = b.lengthSquared();
        if(sq_mag_a < kEpsilon || sq_mag_b < kEpsilon)
        {
            return 0.0f;
        }
        return real_acos(a.dot(b) / (a.length() * b.length()));
    }

    [[nodiscard]] friend TYPE projection(const TYPE& a, const TYPE& b) noexcept
    {
        real sq_mag_b = b.lengthSquared();
        if(sq_mag_b < kEpsilon)
        {
            return TYPE(0.0f);
        }

        /*
         * Basically, vector projection of vector 'a' on a non zero vector 'b' is given by [dot(a, unit(b)) * unit(b)].
         * Checkout https://en.wikipedia.org/wiki/Vector_projection
         */
        return b * (a.dot(b) / sq_mag_b);
    }

    [[nodiscard]] friend TYPE rejection(const TYPE& a, const TYPE& b) noexcept
    {
        /*
         * Rejection of vector 'a' onto vector 'b' is the opposite of projection of vector 'a' onto vector 'b'.
         * To find rejection of 'a' onto 'b', subtract the projection of 'a' onto 'b' from vector 'a'.
         */
        return a - projection(a, b);
    }

    [[nodiscard]] friend TYPE reflection(const TYPE& a, const TYPE& b) noexcept
    {
        return a - (2.0f * projection(a, b));
    }

    [[nodiscard]] friend bool orthonormalize(TYPE& a, TYPE& b, TYPE& c) noexcept
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

    friend std::ostream& operator<<(std::ostream& stream, const TVector& vector)
    {
        stream << "[";
        for(size_type i = 0; i < N; ++i)
        {
            stream << vector.m_data[i];
            if(i < N - 1)
            {
                stream << ", ";
            }
        }
        stream << "]";
        return stream;
    }

protected:
    std::array<T, N> m_data;
};

struct vec2 : public TVector<vec2, real, 2>
{
    real& x = TVector<vec2, real, 2>::m_data[0];
    real& y = TVector<vec2, real, 2>::m_data[1];

    vec2() noexcept            = default;
    vec2(const vec2&) noexcept = default;
    vec2(vec2&&) noexcept      = default;
    vec2(const real& v) noexcept : TVector<vec2, real, 2>(v) {}
    vec2(const real& _x, const real& _y) noexcept : TVector<vec2, real, 2>(_x, _y) {}

    vec2& operator=(const vec2& rhs) noexcept
    {
        if(*this != rhs)
        {
            x = rhs.x;
            y = rhs.y;
        }
        return *this;
    }

    vec2& operator=(vec2&& rhs) noexcept
    {
        if(*this != rhs)
        {
            x = rhs.x;
            y = rhs.y;
        }
        return *this;
    }
};

struct vec3 : public TVector<vec3, real, 3>
{
    real& x = TVector<vec3, real, 3>::m_data[0];
    real& y = TVector<vec3, real, 3>::m_data[1];
    real& z = TVector<vec3, real, 3>::m_data[2];

    vec3() noexcept            = default;
    vec3(const vec3&) noexcept = default;
    vec3(vec3&&) noexcept      = default;
    vec3(const real& v) noexcept : TVector<vec3, real, 3>(v) {}
    vec3(const real& _x, const real& _y, const real& _z) noexcept : TVector<vec3, real, 3>(_x, _y, _z) {}

    vec3& operator=(const vec3& rhs) noexcept
    {
        if(*this != rhs)
        {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
        }
        return *this;
    }

    vec3& operator=(vec3&& rhs) noexcept
    {
        if(*this != rhs)
        {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
        }
        return *this;
    }

    [[nodiscard]] vec2 xy() const noexcept { return {x, y}; }

    [[nodiscard]] vec2 yz() const noexcept { return {y, z}; }

    [[nodiscard]] vec2 xz() const noexcept { return {x, z}; }
};

struct vec4 : public TVector<vec4, real, 4>
{
    real& x = TVector<vec4, real, 4>::m_data[0];
    real& y = TVector<vec4, real, 4>::m_data[1];
    real& z = TVector<vec4, real, 4>::m_data[2];
    real& w = TVector<vec4, real, 4>::m_data[3];

    vec4() noexcept            = default;
    vec4(const vec4&) noexcept = default;
    vec4(vec4&&) noexcept      = default;
    vec4(const real& v) noexcept : TVector<vec4, real, 4>(v) {}
    vec4(const real& _x, const real& _y, const real& _z, const real& _w) noexcept
        : TVector<vec4, real, 4>(_x, _y, _z, _w)
    {
    }

    vec4& operator=(const vec4& rhs) noexcept
    {
        if(*this != rhs)
        {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            w = rhs.w;
        }
        return *this;
    }

    vec4& operator=(vec4&& rhs) noexcept
    {
        if(*this != rhs)
        {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            w = rhs.w;
        }
        return *this;
    }

    [[nodiscard]] vec3 xyz() const noexcept { return {x, y, z}; }

    [[nodiscard]] vec2 xy() const noexcept { return {x, y}; }

    [[nodiscard]] vec2 yz() const noexcept { return {y, z}; }

    [[nodiscard]] vec2 xz() const noexcept { return {x, z}; }
};

} // namespace ramanujan::experimental

#endif