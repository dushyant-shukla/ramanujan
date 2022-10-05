#pragma once

#include "vec3.h"
#include "vec4.h"

namespace ramanujan
{

struct Mat4
{

    bool operator==(const Mat4& other);
    bool operator!=(const Mat4& other);

    friend Mat4 operator+(const Mat4& a, const Mat4& b);
    friend Mat4 operator*(const Mat4& a, float f);
    friend Mat4 operator*(const Mat4& a, const Mat4& b);
    friend Vec4 operator*(const Mat4& m, const Vec4&);

    Vec3 TransformVector(const Mat4& m, const Vec3& v);
    Vec3 TransformPoint(const Mat4& m, const Vec3& v);
    Vec3 TransformPoint(const Mat4& m, const Vec3& v, float& w);

    union
    {
        float v[16];

        // for access to elements based on basis vectors
        struct
        {
            Vec4 right;
            Vec4 up;
            Vec4 forward;
            Vec4 position;
        };

        // for access to elements based on named pairs[basis vector(right, up, forward, position), component(x, y, z,
        // w)]
        struct
        {
            /* col#1*/
            float xx;
            float xy;
            float xz;
            float xw;

            /* col#2*/
            float yx;
            float yy;
            float yz;
            float yw;

            /* col#3*/
            float zx;
            float zy;
            float zz;
            float zw;

            /* col#4*/
            float tx;
            float ty;
            float tz;
            float tw;
        };

        // for element access using column-row notation
        struct
        {

            /* col#1*/
            float c0r0;
            float c0r1;
            float c0r2;
            float c0r3;

            /* col#2*/
            float c1r0;
            float c1r1;
            float c1r2;
            float c1r3;

            /* col#3*/
            float c2r0;
            float c2r1;
            float c2r2;
            float c2r3;

            /* col#4*/
            float c3r0;
            float c3r1;
            float c3r2;
            float c3r3;
        };

        // for element access using row-column notation
        struct
        {

            /* col#1*/
            float r0c0;
            float r1c0;
            float r2c0;
            float r3c0;

            /* col#2*/
            float r0c1;
            float r1c1;
            float r2c1;
            float r3c1;

            /* col#3*/
            float r0c2;
            float r1c2;
            float r2c2;
            float r3c2;

            /* col#4*/
            float r0c3;
            float r1c3;
            float r2c3;
            float r3c3;
        };
    }; // union

    inline Mat4()
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

    inline Mat4(float* fv)
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

    inline Mat4(float _00, // col#1
                float _01,
                float _02,
                float _03,
                float _10, // col#2
                float _11,
                float _12,
                float _13,
                float _20, // col#3
                float _21,
                float _22,
                float _23,
                float _30, // col#4
                float _31,
                float _32,
                float _33)
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

}; // struct Mat4

} // namespace ramanujan