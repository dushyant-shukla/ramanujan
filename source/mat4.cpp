#include "mat4.h"
#include "constants.h"

#include <math.h>

namespace ramanujan
{

bool Mat4::operator==(const Mat4& other)
{
    for(int i = 0; i < 16; ++i)
    {
        if(fabsf(v[i] - other.v[i]) > constants::EPSILON)
        {
            return false;
        }
    }
    return true;
}

bool Mat4::operator!=(const Mat4& other)
{
    return !(*this == other);
}

Mat4 operator+(const Mat4& a, const Mat4& b)
{
    return Mat4(a.xx + b.xx,
                a.xy + b.xy,
                a.xz + b.xz,
                a.xw + b.xw,
                a.yx + b.yx,
                a.yy + b.yy,
                a.yz + b.yz,
                a.yw + b.yw,
                a.zx + b.zx,
                a.zy + b.zy,
                a.zz + b.zz,
                a.zw + b.zw,
                a.tx + b.tx,
                a.ty + b.ty,
                a.tz + b.tz,
                a.tw + b.tw);
}

Mat4 operator*(const Mat4& a, float f)
{
    return Mat4(a.xx * f,
                a.xy * f,
                a.xz * f,
                a.xw * f,
                a.yx * f,
                a.yy * f,
                a.yz * f,
                a.yw * f,
                a.zx * f,
                a.zy * f,
                a.zz * f,
                a.zw * f,
                a.tx * f,
                a.ty * f,
                a.tz * f,
                a.tw * f);
}

Mat4 operator*(const Mat4& a, const Mat4& b)
{
    return Mat4();
}

} // namespace ramanujan