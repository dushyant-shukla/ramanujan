#pragma once

namespace ramanujan
{
template <typename T>
struct TVector2
{
    union
    {
        struct
        {
            T x;
            T y;
        };
        T v[2];
    };

    inline TVector2() : x(T(0)), y(T(0)) {}
    inline TVector2(T v) : x(v), y(v) {}
    inline TVector2(T _x, T _y) : x(_x), y(_y) {}
    inline TVector2(T* fv) : x(fv[0]), y(fv[1]) {}
};

typedef TVector2<float>        Vector2;
typedef TVector2<int>          IVector2;
typedef TVector2<unsigned int> UIVector2;
} // namespace ramanujan
