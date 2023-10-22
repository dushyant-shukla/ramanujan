#pragma once

#include "vector3.h"
#include "vector4.h"

#include "matrix.hpp"

namespace ramanujan
{

struct Matrix4
{
    union
    {
        real v[16];

        // for access to elements based on basis vectors
        struct
        {
            Vector4 right;
            Vector4 up;
            Vector4 forward;
            Vector4 position;
        };

        // for access to elements based on named pairs[basis vector(right, up, forward, position), component(x, y, z,
        // w)]
        struct
        {
            /* col#1*/
            real xx;
            real xy;
            real xz;
            real xw;

            /* col#2*/
            real yx;
            real yy;
            real yz;
            real yw;

            /* col#3*/
            real zx;
            real zy;
            real zz;
            real zw;

            /* col#4*/
            real tx;
            real ty;
            real tz;
            real tw;
        };

        // for element access using column-row notation
        struct
        {

            /* col#1*/
            real c0r0;
            real c0r1;
            real c0r2;
            real c0r3;

            /* col#2*/
            real c1r0;
            real c1r1;
            real c1r2;
            real c1r3;

            /* col#3*/
            real c2r0;
            real c2r1;
            real c2r2;
            real c2r3;

            /* col#4*/
            real c3r0;
            real c3r1;
            real c3r2;
            real c3r3;
        };

        // for element access using row-column notation
        struct
        {

            /* col#1*/
            real r0c0;
            real r1c0;
            real r2c0;
            real r3c0;

            /* col#2*/
            real r0c1;
            real r1c1;
            real r2c1;
            real r3c1;

            /* col#3*/
            real r0c2;
            real r1c2;
            real r2c2;
            real r3c2;

            /* col#4*/
            real r0c3;
            real r1c3;
            real r2c3;
            real r3c3;
        };
    }; // union

    inline Matrix4()
        : xx(1) // col#1
        , xy(0)
        , xz(0)
        , xw(0)
        , yx(0) // col#2
        , yy(1)
        , yz(0)
        , yw(0)
        , zx(0) // col#3
        , zy(0)
        , zz(1)
        , zw(0)
        , tx(0) // col#4
        , ty(0)
        , tz(0)
        , tw(1)
    {
    }

    inline Matrix4(real* fv)
        : xx(fv[0]) // col#1
        , xy(fv[1])
        , xz(fv[2])
        , xw(fv[3])
        , yx(fv[4]) // col#2
        , yy(fv[5])
        , yz(fv[6])
        , yw(fv[7])
        , zx(fv[8]) // col#3
        , zy(fv[9])
        , zz(fv[10])
        , zw(fv[11])
        , tx(fv[12]) // col#4
        , ty(fv[13])
        , tz(fv[14])
        , tw(fv[15])
    {
    }

    inline Matrix4(real _00, // col#1
                   real _01,
                   real _02,
                   real _03,
                   real _10, // col#2
                   real _11,
                   real _12,
                   real _13,
                   real _20, // col#3
                   real _21,
                   real _22,
                   real _23,
                   real _30, // col#4
                   real _31,
                   real _32,
                   real _33)
        : xx(_00)
        , xy(_01)
        , xz(_02)
        , xw(_03)
        , yx(_10)
        , yy(_11)
        , yz(_12)
        , yw(_13)
        , zx(_20)
        , zy(_21)
        , zz(_22)
        , zw(_23)
        , tx(_30)
        , ty(_31)
        , tz(_32)
        , tw(_33)
    {
    }

    bool           operator==(const Matrix4& other);
    bool           operator!=(const Matrix4& other);
    friend Matrix4 operator+(const Matrix4& a, const Matrix4& b);
    friend Matrix4 operator*(const Matrix4& a, real f);
    friend Matrix4 operator*(const Matrix4& a, const Matrix4& b);
    friend Vector4 operator*(const Matrix4& m, const Vector4&);

    operator ramanujan::experimental::mat4() const
    {
        return ramanujan::experimental::mat4{xx, xy, xz, xw, yx, yy, yz, yw, zx, zy, zz, zw, tx, ty, tz, tw};
    }

}; // struct Mat4

void    Transpose(Matrix4& m);
Matrix4 Transposed(const Matrix4& m);
Vector3 TransformVector(const Matrix4& m, const Vector3& v);
Vector3 TransformPoint(const Matrix4& m, const Vector3& v);
Vector3 TransformPoint(const Matrix4& m, const Vector3& v, real& w);
real    Determinant(const Matrix4& m);
Matrix4 Adjugate(const Matrix4& m);
Matrix4 Inverse(const Matrix4& m);
void    Invert(Matrix4& m);

} // namespace ramanujan