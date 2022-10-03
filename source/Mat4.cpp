#include "Mat4.h"
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
    return Mat4();
}

Mat4 operator*(const Mat4& a, float f)
{
    return Mat4();
}

Mat4 operator*(const Mat4& a, const Mat4& b)
{
    return Mat4();
}

} // namespace ramanujan