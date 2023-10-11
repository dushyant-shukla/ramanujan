#ifndef RAMANUJAN_VECTOR
#define RAMANUJAN_VECTOR

#include "constants.h"
#include "precision.h"

#include <array>

namespace ramanujan::experimental
{

template <typename T, unsigned int DIMENSION>
class Vector
{

public:
    Vector();

    Vector(const T& value);

    template <typename... Args>
    Vector(Args... args);

    T x() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 2, T>::type y() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 3, T>::type z() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 4, T>::type w() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 4, Vector<T, 3>>::type xyz() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 3, Vector<T, 2>>::type xy() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 3, Vector<T, 2>>::type xz() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 3, Vector<T, 2>>::type yz() const;

    T r() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 2, T>::type g() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 3, T>::type b() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 4, T>::type a() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 4, Vector<T, 3>>::type rgb() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 3, Vector<T, 2>>::type rg() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 3, Vector<T, 2>>::type rb() const;

    template <unsigned int D = DIMENSION>
    typename std::enable_if<D >= 3, Vector<T, 2>>::type gb() const;

    T const* data() const;
    T*       data();

    T&       operator[](unsigned int index);
    T const& operator[](unsigned int index) const;
    bool     operator==(T const& rhs) const;
    bool     operator!=(T const& rhs) const;

    Vector<T, DIMENSION>& operator++();
    Vector<T, DIMENSION>  operator++(int);

    Vector<T, DIMENSION>& operator--();
    Vector<T, DIMENSION>  operator--(int);

    Vector<T, DIMENSION>& operator+=(Vector<T, DIMENSION> const& rhs);

    template <typename U>
    Vector<T, DIMENSION>& operator+=(U scalar);

    Vector<T, DIMENSION>& operator-=(Vector<T, DIMENSION> const& rhs);

    template <typename U>
    Vector<T, DIMENSION>& operator-=(U scalar);

    Vector<T, DIMENSION>& operator*=(Vector<T, DIMENSION> const& rhs);

    template <typename U>
    Vector<T, DIMENSION>& operator*=(U scalar);

    Vector<T, DIMENSION>& operator/=(Vector<T, DIMENSION> const& rhs);

    template <typename U>
    Vector<T, DIMENSION>& operator/=(U scalar);

    friend Vector<T, DIMENSION> operator+(Vector<T, DIMENSION> const& lhs, Vector<T, DIMENSION> const& rhs);
    friend Vector<T, DIMENSION> operator+(Vector<T, DIMENSION> const& lhs, T scalar);
    friend Vector<T, DIMENSION> operator+(T scalar, Vector<T, DIMENSION> const& rhs);

    friend Vector<T, DIMENSION> operator-(Vector<T, DIMENSION> const& lhs, Vector<T, DIMENSION> const& rhs);
    friend Vector<T, DIMENSION> operator-(Vector<T, DIMENSION> const& lhs, T scalar);
    friend Vector<T, DIMENSION> operator-(T scalar, Vector<T, DIMENSION> const& rhs);

    friend Vector<T, DIMENSION> operator*(Vector<T, DIMENSION> const& lhs, Vector<T, DIMENSION> const& rhs);
    friend Vector<T, DIMENSION> operator*(Vector<T, DIMENSION> const& lhs, T scalar);
    friend Vector<T, DIMENSION> operator*(T scalar, Vector<T, DIMENSION> const& rhs);

    friend Vector<T, DIMENSION> operator/(Vector<T, DIMENSION> const& lhs, Vector<T, DIMENSION> const& rhs);
    friend Vector<T, DIMENSION> operator/(Vector<T, DIMENSION> const& lhs, T scalar);
    friend Vector<T, DIMENSION> operator/(T scalar, Vector<T, DIMENSION> const& rhs);

    real                 dot(Vector<T, DIMENSION> const& rhs) const;
    Vector<T, DIMENSION> cross(Vector<T, DIMENSION> const& rhs) const;
    Vector<T, DIMENSION> reflect(Vector<T, DIMENSION> const& rhs) const;
    Vector<T, DIMENSION> reject(Vector<T, DIMENSION> const& rhs) const;
    Vector<T, DIMENSION> angle(Vector<T, DIMENSION> const& rhs) const;
    Vector<T, DIMENSION> project(Vector<T, DIMENSION> const& rhs) const;
    real                 length() const;
    real                 lengthSquare() const;
    real                 distance(Vector<T, DIMENSION> const& rhs) const;
    void                 normalize();
    Vector<T, DIMENSION> normalize() const;

    friend Vector<T, DIMENSION> lerp(Vector<T, DIMENSION> const& start, Vector<T, DIMENSION> const& end, real t);
    friend Vector<T, DIMENSION> nlerp(Vector<T, DIMENSION> const& start, Vector<T, DIMENSION> const& end, real t);
    friend Vector<T, DIMENSION> slerp(Vector<T, DIMENSION> const& start, Vector<T, DIMENSION> const& end, real t);
    friend bool makeOrthonormalBasis(Vector<T, DIMENSION>& a, Vector<T, DIMENSION>& b, Vector<T, DIMENSION>& c);

private:
    std::array<T, DIMENSION> m_data;
};

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION>::Vector()
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        m_data[i] = T(0);
    }
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION>::Vector(const T& value)
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        m_data[i] = value;
    }
}

template <typename T, unsigned int DIMENSION>
template <typename... Args>
inline Vector<T, DIMENSION>::Vector(Args... args) : m_data{args...}
{
}

template <typename T, unsigned int DIMENSION>
inline T Vector<T, DIMENSION>::x() const
{
    return m_data[0];
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 2, T>::type Vector<T, DIMENSION>::y() const
{
    return m_data[1];
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 3, T>::type Vector<T, DIMENSION>::z() const
{
    return m_data[2];
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 4, T>::type Vector<T, DIMENSION>::w() const
{
    return m_data[3];
}

template <typename T, unsigned int DIMENSION>
template <unsigned D>
inline typename std::enable_if<D >= 4, Vector<T, 3>>::type Vector<T, DIMENSION>::xyz() const
{
    return Vector<T, 3>{m_data[0], m_data[1], m_data[2]};
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 3, Vector<T, 2>>::type Vector<T, DIMENSION>::xy() const
{
    return Vector<T, 2>{m_data[0], m_data[1]};
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 3, Vector<T, 2>>::type Vector<T, DIMENSION>::xz() const
{
    return Vector<T, 2>{m_data[0], m_data[2]};
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 3, Vector<T, 2>>::type Vector<T, DIMENSION>::yz() const
{
    return Vector<T, 2>{m_data[1], m_data[2]};
}

template <typename T, unsigned int DIMENSION>
inline T Vector<T, DIMENSION>::r() const
{
    return m_data[0];
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 2, T>::type Vector<T, DIMENSION>::g() const
{
    return m_data[1];
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 3, T>::type Vector<T, DIMENSION>::b() const
{
    return m_data[2];
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 4, T>::type Vector<T, DIMENSION>::a() const
{
    return m_data[3];
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 4, Vector<T, 3>>::type Vector<T, DIMENSION>::rgb() const
{
    return Vector<T, 3>{m_data[0], m_data[1], m_data[2]};
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 3, Vector<T, 2>>::type Vector<T, DIMENSION>::rg() const
{
    return Vector<T, 2>{m_data[0], m_data[1]};
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 3, Vector<T, 2>>::type Vector<T, DIMENSION>::rb() const
{
    return Vector<T, 2>{m_data[0], m_data[2]};
}

template <typename T, unsigned int DIMENSION>
template <unsigned int D>
inline typename std::enable_if<D >= 3, Vector<T, 2>>::type Vector<T, DIMENSION>::gb() const
{
    return Vector<T, 2>{m_data[1], m_data[2]};
}

template <typename T, unsigned int DIMENSION>
inline T const* Vector<T, DIMENSION>::data() const
{
    return m_data.data();
}

template <typename T, unsigned int DIMENSION>
inline T* Vector<T, DIMENSION>::data()
{
    return m_data.data();
}

template <typename T, unsigned int DIMENSION>
inline T& Vector<T, DIMENSION>::operator[](unsigned int index)
{
    return m_data[index];
}

template <typename T, unsigned int DIMENSION>
inline T const& Vector<T, DIMENSION>::operator[](unsigned int index) const
{
    return m_data[index];
}

template <typename T, unsigned int DIMENSION>
inline bool Vector<T, DIMENSION>::operator==(T const& rhs) const
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        if(this->m_data[i] != rhs.m_data[i])
        {
            return false;
        }
    }
    return true;
}

template <typename T, unsigned int DIMENSION>
inline bool Vector<T, DIMENSION>::operator!=(T const& rhs) const
{
    return !(*this == rhs);
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION>& Vector<T, DIMENSION>::operator++()
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        ++(this->m_data[i]);
    }
    return *this;
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION> Vector<T, DIMENSION>::operator++(int)
{
    Vector<T, DIMENSION> result(*this);
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        ++result.m_data[i];
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION>& Vector<T, DIMENSION>::operator--()
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        ++(this->m_data[i]);
    }
    return *this;
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION> Vector<T, DIMENSION>::operator--(int)
{
    Vector<T, DIMENSION> result(*this);
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        --result.m_data[i];
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION>& Vector<T, DIMENSION>::operator+=(Vector<T, DIMENSION> const& rhs)
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        this->m_data[i] += rhs.m_data[i];
    }
    return *this;
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION>& Vector<T, DIMENSION>::operator-=(Vector<T, DIMENSION> const& rhs)
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        this->m_data[i] -= rhs.m_data[i];
    }
    return *this;
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION>& Vector<T, DIMENSION>::operator*=(Vector<T, DIMENSION> const& rhs)
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        this->m_data[i] *= rhs.m_data[i];
    }
    return *this;
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION>& Vector<T, DIMENSION>::operator/=(Vector<T, DIMENSION> const& rhs)
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        this->m_data[i] /= rhs.m_data[i];
    }
    return *this;
}

template <typename T, unsigned int DIMENSION>
template <typename U>
inline Vector<T, DIMENSION>& Vector<T, DIMENSION>::operator+=(U scalar)
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        this->m_data[i] += static_cast<T>(scalar);
    }
    return *this;
}

template <typename T, unsigned int DIMENSION>
template <typename U>
inline Vector<T, DIMENSION>& Vector<T, DIMENSION>::operator-=(U scalar)
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        this->m_data[i] -= static_cast<T>(scalar);
    }
    return *this;
}

template <typename T, unsigned int DIMENSION>
template <typename U>
inline Vector<T, DIMENSION>& Vector<T, DIMENSION>::operator*=(U scalar)
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        this->m_data[i] *= static_cast<T>(scalar);
    }
    return *this;
}

template <typename T, unsigned int DIMENSION>
template <typename U>
inline Vector<T, DIMENSION>& Vector<T, DIMENSION>::operator/=(U scalar)
{
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        this->m_data[i] /= static_cast<T>(scalar);
    }
    return *this;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> operator+(Vector<T, DIMENSION> const& lhs, Vector<T, DIMENSION> const& rhs)
{
    Vector<T, DIMENSION> result;
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        result.m_data[i] = static_cast<T>(lhs.m_data[i] + rhs.m_data[i]);
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> operator+(Vector<T, DIMENSION> const& lhs, T scalar)
{
    Vector<T, DIMENSION> result;
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        result.m_data[i] = lhs.m_data[i] + scalar;
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> operator+(T scalar, Vector<T, DIMENSION> const& rhs)
{
    return rhs + scalar;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> operator-(Vector<T, DIMENSION> const& lhs, Vector<T, DIMENSION> const& rhs)
{
    Vector<T, DIMENSION> result;
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        result.m_data[i] = lhs.m_data[i] - rhs.m_data[i];
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> operator-(Vector<T, DIMENSION> const& lhs, T scalar)
{
    Vector<T, DIMENSION> result;
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        result.m_data[i] = lhs.m_data[i] - scalar;
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> operator-(T scalar, Vector<T, DIMENSION> const& rhs)
{
    return rhs - scalar;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> operator*(Vector<T, DIMENSION> const& lhs, Vector<T, DIMENSION> const& rhs)
{
    Vector<T, DIMENSION> result;
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        result.m_data[i] = static_cast<T>(lhs.m_data[i] * rhs.m_data[i]);
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> operator*(Vector<T, DIMENSION> const& lhs, T scalar)
{
    Vector<T, DIMENSION> result;
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        result.m_data[i] = static_cast<T>(scalar * rhs.m_data[i]);
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> operator*(T scalar, Vector<T, DIMENSION> const& rhs)
{
    return rhs * scalar;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> operator/(Vector<T, DIMENSION> const& lhs, Vector<T, DIMENSION> const& rhs)
{
    Vector<T, DIMENSION> result;
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        result.m_data[i] = static_cast<T>(lhs.m_data[i] / rhs.m_data[i]);
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> operator/(Vector<T, DIMENSION> const& lhs, T scalar)
{
    Vector<T, DIMENSION> result;
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        result.m_data[i] = static_cast<T>(lhs.m_data[i] / scalar);
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> operator/(T scalar, Vector<T, DIMENSION> const& rhs)
{
    return rhs / scalar;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> lerp(Vector<T, DIMENSION> const& start, Vector<T, DIMENSION> const& end, real t)
{
    Vector<T, DIMENSION> result;
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        result[i] = start[i] + (end[i] - start[i]) * t;
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> nlerp(Vector<T, DIMENSION> const& start, Vector<T, DIMENSION> const& end, real t)
{
    Vector<T, DIMENSION> result = lerp(start, end, t);
    return result.normalize();
}

template <typename T, unsigned int DIMENSION>
Vector<T, DIMENSION> slerp(Vector<T, DIMENSION> const& start, Vector<T, DIMENSION> const& end, real t)
{
    if(t < 0.01f)
    {
        return lerp(start, end, t);
    }

    auto from = start.normalize();
    auto to   = end.normalize();

    real theta     = from.angle(to);
    real sin_theta = sinf(theta);

    real a = sinf((1.0f - t) * theta) / sin_theta;
    real b = sinf(t * theta) / sin_theta;

    return a * from + b * to;
}

template <typename T, unsigned int DIMENSION>
bool makeOrthonormalBasis(Vector<T, DIMENSION>& a, Vector<T, DIMENSION>& b, Vector<T, DIMENSION>& c)
{
    // TODO: Should this be part of the math library?
    /*
     * Note that the construction of an orthonormal basis is a situation where it matters a great deal whether you are
     * working in a left - or right - handed coordinate system. The following algorithm is designed for right - handed
     * systems. If you need a left - handed coordinate system, then you can simply change the order of the operands for
     * both the cross - products. This will give you a left - handed orthonormal basis.
     */
    a.normalize();
    c             = a.cross(b); // change to Cross(b, a) for a left-handed coordinate systems
    real c_mag_sq = c.lengthSquare();
    if(c_mag_sq < Constants::EPSILON)
    {
        return false;
    }
    c.normalize();
    b = c.cross(a); // change to Cross(a, c) for a left-handed coordinate systems
    return true;
}

using Vec2 = Vector<float, 2>;
using Vec3 = Vector<float, 3>;
using Vec4 = Vector<float, 4>;

using Color3 = Vector<float, 3>;
using Color4 = Vector<float, 4>;

using IVec2 = Vector<int, 2>;
using IVec3 = Vector<int, 3>;
using IVec4 = Vector<int, 4>;

using UVec2 = Vector<unsigned int, 2>;
using UVec3 = Vector<unsigned int, 3>;
using UVec4 = Vector<unsigned int, 4>;

template <typename T, unsigned int DIMENSION>
inline real Vector<T, DIMENSION>::dot(Vector<T, DIMENSION> const& rhs) const
{
    real result = 0.0f;
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        result += static_cast<T>(lhs[i] * rhs[i]);
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION> Vector<T, DIMENSION>::cross(Vector<T, DIMENSION> const& rhs) const
{
    return Vector<T, DIMENSION>();
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION> Vector<T, DIMENSION>::reflect(Vector<T, DIMENSION> const& rhs) const
{
    return Vector<T, DIMENSION>();
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION> Vector<T, DIMENSION>::reject(Vector<T, DIMENSION> const& rhs) const
{
    return Vector<T, DIMENSION>();
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION> Vector<T, DIMENSION>::angle(Vector<T, DIMENSION> const& rhs) const
{
    return Vector<T, DIMENSION>();
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION> Vector<T, DIMENSION>::project(Vector<T, DIMENSION> const& rhs) const
{
    return Vector<T, DIMENSION>();
}

template <typename T, unsigned int DIMENSION>
inline real Vector<T, DIMENSION>::length() const
{
    real length_square = lengthSquare();
    if(length_sq < Constants::EPSILON)
    {
        return 0.0f;
    }
    return sqrtf(length_sq);
}

template <typename T, unsigned int DIMENSION>
inline real Vector<T, DIMENSION>::lengthSquare() const
{
    real result = 0.0f;
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        result += static_cast<T>(*this[i] * *this[i]);
    }
    return result;
}

template <typename T, unsigned int DIMENSION>
inline real Vector<T, DIMENSION>::distance(Vector<T, DIMENSION> const& rhs) const
{
    return this->length(rhs);
}

template <typename T, unsigned int DIMENSION>
inline void Vector<T, DIMENSION>::normalize()
{
    real length_sq = lengthSquare();
    if(length_sq < Constants::EPSILON)
    {
        return;
    }
    real inverted_length = 1.0f / sqrtf(length_sq);
    for(unsigned int i = 0; i < DIMENSION; ++i)
    {
        this->m_data[i] *= inverted_length;
    }
}

template <typename T, unsigned int DIMENSION>
inline Vector<T, DIMENSION> Vector<T, DIMENSION>::normalize() const
{
    Vector<T, DIMENSION> result;
    real                 length_sq = lengthSquare();
    if(length_sq < Constants::EPSILON)
    {
        return *this;
    }
    real inverted_length = 1.0f / sqrtf(length_sq);
    result *= inverted_length;
    // for(unsigned int i = 0; i < DIMENSION; ++i)
    //{
    //     result[i] *= inverted_length;
    // }
    return result;
}

} // namespace ramanujan::experimental

#endif
