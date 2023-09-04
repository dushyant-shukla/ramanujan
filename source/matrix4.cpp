#include "matrix4.h"
#include "constants.h"

#include <math.h>

namespace ramanujan
{

bool Matrix4::operator==(const Matrix4& other)
{
    for(int i = 0; i < 16; ++i)
    {
        if(fabsf(v[i] - other.v[i]) > Constants::EPSILON)
        {
            return false;
        }
    }
    return true;
}

bool Matrix4::operator!=(const Matrix4& other)
{
    return !(*this == other);
}

Matrix4 operator+(const Matrix4& a, const Matrix4& b)
{
    return Matrix4(a.xx + b.xx,
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

Matrix4 operator*(const Matrix4& a, real f)
{
    return Matrix4(a.xx * f,
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

Matrix4 operator*(const Matrix4& a, const Matrix4& b)
{
    return Matrix4(M4_DOT(0, 0),
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
#define M4_V4_DOT(row_index, x, y, z, w) \
    m.v[0 * 4 + row_index] * x + m.v[1 * 4 + row_index] * y + m.v[2 * 4 + row_index] * z + m.v[3 * 4 + row_index] * w

Vector4 operator*(const Matrix4& m, const Vector4& v)
{
    return Vector4(M4_V4_DOT(0, v.x, v.y, v.z, v.w),
                   M4_V4_DOT(1, v.x, v.y, v.z, v.w),
                   M4_V4_DOT(2, v.x, v.y, v.z, v.w),
                   M4_V4_DOT(3, v.x, v.y, v.z, v.w));
}

#define M4_SWAP(x, y) \
    {                 \
        real t = x;   \
        x      = y;   \
        y      = t;   \
    }

/**
 * Transposing a matrix is useful if you need to convert a matrix from row-major to column-major or vice versa.
 *
 * \param m Matrix to be transposed
 */
void Transpose(Matrix4& m)
{
    M4_SWAP(m.yx, m.xy);
    M4_SWAP(m.zx, m.xz);
    M4_SWAP(m.tx, m.xw);
    M4_SWAP(m.zy, m.yz);
    M4_SWAP(m.ty, m.yw);
    M4_SWAP(m.tz, m.zw);
}

/**
 * Transposing a matrix is useful if you need to convert a matrix from row-major to column-major or vice versa.
 *
 * \param m Matrix to be transposed
 * \return Transposed matrix
 */
Matrix4 Transposed(const Matrix4& m)
{
    // Transpose is converting a column-major matrix to a row-major matrix here
    return Matrix4(m.xx, m.yx, m.zx, m.tx, m.xy, m.yy, m.zy, m.ty, m.xz, m.yz, m.zz, m.tz, m.xw, m.yw, m.zw, m.tw);
}

Vector3 TransformVector(const Matrix4& m, const Vector3& v)
{
    return Vector3(M4_V4_DOT(0, v.x, v.y, v.z, 0.0f),
                   M4_V4_DOT(1, v.x, v.y, v.z, 0.0f),
                   M4_V4_DOT(2, v.x, v.y, v.z, 0.0f));
}

Vector3 TransformPoint(const Matrix4& m, const Vector3& v)
{
    return Vector3(M4_V4_DOT(0, v.x, v.y, v.z, 1.0f),
                   M4_V4_DOT(1, v.x, v.y, v.z, 1.0f),
                   M4_V4_DOT(2, v.x, v.y, v.z, 1.0f));
}

Vector3 TransformPoint(const Matrix4& m, const Vector3& v, real& w)
{
    real _w = w;
    w       = M4_V4_DOT(3, v.x, v.y, v.z, _w);
    return Vector3(M4_V4_DOT(0, v.x, v.y, v.z, _w), M4_V4_DOT(1, v.x, v.y, v.z, _w), M4_V4_DOT(2, v.x, v.y, v.z, _w));
}

#define M4_3X3MINOR(x, c0, c1, c2, r0, r1, r2)                                              \
    (x[c0 * 4 + r0] * (x[c1 * 4 + r1] * x[c2 * 4 + r2] - x[c1 * 4 + r2] * x[c2 * 4 + r1]) - \
     x[c1 * 4 + r0] * (x[c0 * 4 + r1] * x[c2 * 4 + r2] - x[c0 * 4 + r2] * x[c2 * 4 + r1]) + \
     x[c2 * 4 + r0] * (x[c0 * 4 + r1] * x[c1 * 4 + r2] - x[c0 * 4 + r2] * x[c1 * 4 + r1]))

real Determinant(const Matrix4& m)
{
    return m.v[0] * M4_3X3MINOR(m.v, 1, 2, 3, 1, 2, 3) - m.v[4] * M4_3X3MINOR(m.v, 0, 2, 3, 1, 2, 3) +
           m.v[8] * M4_3X3MINOR(m.v, 0, 1, 3, 1, 2, 3) - m.v[12] * M4_3X3MINOR(m.v, 0, 1, 2, 1, 2, 3);
}

Matrix4 Adjugate(const Matrix4& m)
{
    // Cof (M[i, j]) = Minor(M[i, j]] * pow(-1, i + j)
    // Having C1, C2, C3, R1, R2, R3 etc constants would be good instead of these raw indices
    Matrix4 cofactor;
    cofactor.v[0]  = M4_3X3MINOR(m.v, 1, 2, 3, 1, 2, 3);
    cofactor.v[1]  = -M4_3X3MINOR(m.v, 1, 2, 3, 0, 2, 3);
    cofactor.v[2]  = M4_3X3MINOR(m.v, 1, 2, 3, 0, 1, 3);
    cofactor.v[3]  = -M4_3X3MINOR(m.v, 1, 2, 3, 0, 1, 2);
    cofactor.v[4]  = -M4_3X3MINOR(m.v, 0, 2, 3, 1, 2, 3);
    cofactor.v[5]  = M4_3X3MINOR(m.v, 0, 2, 3, 0, 2, 3);
    cofactor.v[6]  = -M4_3X3MINOR(m.v, 0, 2, 3, 0, 1, 3);
    cofactor.v[7]  = M4_3X3MINOR(m.v, 0, 2, 3, 0, 1, 2);
    cofactor.v[8]  = M4_3X3MINOR(m.v, 0, 1, 3, 1, 2, 3);
    cofactor.v[9]  = -M4_3X3MINOR(m.v, 0, 1, 3, 0, 2, 3);
    cofactor.v[10] = M4_3X3MINOR(m.v, 0, 1, 3, 0, 1, 3);
    cofactor.v[11] = -M4_3X3MINOR(m.v, 0, 1, 3, 0, 1, 2);
    cofactor.v[12] = -M4_3X3MINOR(m.v, 0, 1, 2, 1, 2, 3);
    cofactor.v[13] = M4_3X3MINOR(m.v, 0, 1, 2, 0, 2, 3);
    cofactor.v[14] = -M4_3X3MINOR(m.v, 0, 1, 2, 0, 1, 3);
    cofactor.v[15] = M4_3X3MINOR(m.v, 0, 1, 2, 0, 1, 2);
    return Transposed(cofactor);
}

Matrix4 Inverse(const Matrix4& m)
{
    real det = Determinant(m);
    if(det == 0.0f)
    {
        return Matrix4();
    }
    Matrix4 adj = Adjugate(m);
    return adj * (1.0f / det);
}

void Invert(Matrix4& m)
{
    real det = Determinant(m);
    if(det == 0.0f)
    {
        m = Matrix4();
        return;
    }
    m = Adjugate(m) * (1.0f / det);
}

} // namespace ramanujan