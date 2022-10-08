#pragma once

namespace ramanujan
{
template <typename T>
struct TVector4
{
    union
    {
        struct
        {
            T x;
            T y;
            T z;
            T w;
        };
        T v[4];
    };

    inline TVector4() : x(T(0)), y(T(0)), z(T(0)), w(T(0)) {}
    inline TVector4(T v) : x(v), y(v), z(v), w(v) {}
    inline TVector4(T _x, T _y, T _z, T _w) : x(_x), y(_y), z(_z), w(_w) {}
    inline TVector4(T* fv) : x(fv[0]), y(fv[1]), z(fv[2]), w(fv[3]) {}
};

typedef TVector4<float>        Vector4;
typedef TVector4<int>          IVector4;
typedef TVector4<unsigned int> UIVector4;
} // namespace ramanujan