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
        for(size_type i = 0; i < size(); ++i)
        {
            m_data[i] += rhs.m_data[i];
        }
        return *this;
    }

    constexpr TVector& operator-=(const TVector& rhs) noexcept
    {
        for(size_type i = 0; i < size(); ++i)
        {
            m_data[i] -= rhs.m_data[i];
        }
        return *this;
    }

    constexpr TVector& operator*=(const TVector& rhs) noexcept
    {
        for(size_type i = 0; i < size(); ++i)
        {
            m_data[i] *= rhs.m_data[i];
        }
        return *this;
    }

    constexpr TVector& operator/=(const TVector& rhs) noexcept
    {
        for(size_type i = 0; i < size(); ++i)
        {
            m_data[i] /= rhs.m_data[i];
        }
        return *this;
    }

    constexpr TVector& operator+=(const T& rhs) noexcept
    {
        for(size_type i = 0; i < size(); ++i)
        {
            m_data[i] += rhs;
        }
        return *this;
    }

    constexpr TVector& operator-=(const T& rhs) noexcept
    {
        for(size_type i = 0; i < size(); ++i)
        {
            m_data[i] -= rhs;
        }
        return *this;
    }

    constexpr TVector& operator*=(const T& rhs) noexcept
    {
        for(size_type i = 0; i < size(); ++i)
        {
            m_data[i] *= rhs;
        }
        return *this;
    }

    constexpr TVector& operator/=(const T& rhs) noexcept
    {
        for(size_type i = 0; i < size(); ++i)
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

using color3 = TVec3<float>;
using color4 = TVec4<float>;

} // namespace ramanujan::experimental

#endif