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

// In column-major matrix, c[r, c] = c * R + r
// Each element in the matrix multiplication can be described as the dot product between the corresponding row of the
// first matrix, and the corresponding column of the second matrix. In  short, M[i, j] = Ri(A) dot Ci(B)
// Here we access the elements of the given row from the first matrix, and the elements from the given column from the
// second matrix
#define M4_DOT(a_row, b_col)                                                            \
    a.v[0 * 4 + a_row] * b.v[b_col * 4 + 0] + a.v[1 * 4 + a_row] * b.v[b_col * 4 + 1] + \
        a.v[2 * 4 + a_row] * b.v[b_col * 4 + 2] + a.v[3 * 4 + a_row] * b.v[b_col * 4 + 3]

Mat4 operator*(const Mat4& a, const Mat4& b)
{
    return Mat4(M4_DOT(0, 0),
                M4_DOT(1, 0),
                M4_DOT(2, 0),
                M4_DOT(3, 0),
                M4_DOT(0, 1),
                M4_DOT(1, 1),
                M4_DOT(2, 1),
                M4_DOT(3, 1),
                M4_DOT(0, 2),
                M4_DOT(1, 2),
                M4_DOT(2, 2),
                M4_DOT(3, 2),
                M4_DOT(0, 3),
                M4_DOT(1, 3),
                M4_DOT(2, 3),
                M4_DOT(3, 3));
}

// The 4D vector being transformed can be thought of as a matrix with 4 columns and 1 row.
// When a matrix transforms a vector, it affects both the orientation and scale of the vector (w component of a vector
// is 0). When a matrix transforms a point, it just translates the point in space (w component of a point is 1).
#define M4_V4_DOT(m_row, x, y, z, w) \
    m.v[0 * 4 + m_row] * x + m.v[1 * 4 + m_row] * y + m.v[2 * 4 + m_row] * z + m.v[3 * 4 + m_row] * w

Vec4 operator*(const Mat4& m, const Vec4& v)
{
    return Vec4(M4_V4_DOT(0, v.x, v.y, v.z, v.w),
                M4_V4_DOT(1, v.x, v.y, v.z, v.w),
                M4_V4_DOT(2, v.x, v.y, v.z, v.w),
                M4_V4_DOT(3, v.x, v.y, v.z, v.w));
}

Vec3 Mat4::TransformVector(const Mat4& m, const Vec3& v)
{
    return Vec3(M4_V4_DOT(0, v.x, v.y, v.z, 0.0f),
                M4_V4_DOT(1, v.x, v.y, v.z, 0.0f),
                M4_V4_DOT(2, v.x, v.y, v.z, 0.0f));
}

Vec3 Mat4::TransformPoint(const Mat4& m, const Vec3& v)
{
    return Vec3(M4_V4_DOT(0, v.x, v.y, v.z, 1.0f),
                M4_V4_DOT(1, v.x, v.y, v.z, 1.0f),
                M4_V4_DOT(2, v.x, v.y, v.z, 1.0f));
}

Vec3 Mat4::TransformPoint(const Mat4& m, const Vec3& v, float& w)
{
    float _w = w;
    w        = M4_V4_DOT(3, v.x, v.y, v.z, _w);
    return Vec3(M4_V4_DOT(0, v.x, v.y, v.z, _w), M4_V4_DOT(1, v.x, v.y, v.z, _w), M4_V4_DOT(2, v.x, v.y, v.z, _w));
}

} // namespace ramanujan