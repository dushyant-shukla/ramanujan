#ifndef RAMANUJAN_MATRIX
#define RAMANUJAN_MATRIX

#include "constants.h"
#include "precision.h"
#include "vector.hpp"

#include <array>

namespace ramanujan::experimental
{

template <typename MAT_TYPE, typename T, std::size_t ROWS, std::size_t COLUMNS>
class TMatrix
{
public:
    using value_type      = T;
    using size_type       = std::size_t;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using pointer         = value_type*;
    using const_pointer   = const value_type*;

    static constexpr size_t size = ROWS * COLUMNS;

    MAT_TYPE&       type() { return static_cast<MAT_TYPE&>(*this); }
    const MAT_TYPE& type() const { return static_cast<const MAT_TYPE&>(*this); }

    // Todo:: Not sure is this would work (statics are not inherited)
    static MAT_TYPE identity() noexcept { return MAT_TYPE{}; }

    constexpr size_type rows() const noexcept { return ROWS; }

    constexpr size_type columns() const noexcept { return COLUMNS; }

    // constexpr size_type size() const noexcept { return ROWS * COLUMNS; }

    constexpr pointer data() const noexcept { return type().m.data(); }

    constexpr const_pointer data() noexcept { return type().m.data(); }

    constexpr reference operator()(size_type row, size_type column) noexcept { return type().m[column * ROWS + row]; }

    constexpr const_reference operator()(size_type row, size_type column) const noexcept
    {
        return type().m[column * ROWS + row]; // m[stride + offset]
    }

    MAT_TYPE& operator+=(const MAT_TYPE& rhs) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < size; ++i)
        {
            self.m[i] += rhs.m[i];
        }
        return self;
    }

    MAT_TYPE& operator-=(const MAT_TYPE& rhs) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < size; ++i)
        {
            self.m[i] -= rhs.m[i];
        }
        return self;
    }

    MAT_TYPE& operator/=(const MAT_TYPE& rhs) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < size; ++i)
        {
            assert(rhs.m[i] > T(0));
            self.m[i] /= rhs.m[i];
        }
        return self;
    }

    MAT_TYPE& operator+=(const_reference scalar) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < size; ++i)
        {
            self.m[i] += scalar;
        }
        return self;
    }

    MAT_TYPE& operator-=(const_reference scalar) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < size; ++i)
        {
            self.m[i] -= scalar;
        }
        return self;
    }

    MAT_TYPE& operator*=(const_reference scalar) noexcept
    {
        auto& self = type();
        for(size_type i = 0; i < size; ++i)
        {
            self.m[i] *= scalar;
        }
        return self;
    }

    MAT_TYPE& operator/=(const_reference scalar) noexcept
    {
        assert(scalar > T(0));
        auto& self = type();
        for(size_type i = 0; i < size; ++i)
        {
            self.m[i] /= scalar;
        }
        return self;
    }

    bool operator==(const MAT_TYPE& rhs) const noexcept
    {
        const auto& self = type();
        for(size_type i = 0; i < size; ++i)
        {
            if(real_abs(self.m[i] - rhs.m[i]) > T(kEpsilon))
            {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const MAT_TYPE& rhs) const noexcept { return !(type() == rhs); }

    MAT_TYPE operator+(const MAT_TYPE& rhs) const noexcept
    {
        auto&    self = type();
        MAT_TYPE result{};
        for(size_type i = 0; i < size; ++i)
        {
            result.m[i] = self.m[i] + rhs.m[i];
        }
        return result;
    }

    MAT_TYPE operator-(const MAT_TYPE& rhs) const noexcept
    {
        auto&    self = type();
        MAT_TYPE result{};
        for(size_type i = 0; i < size; ++i)
        {
            result.m[i] = self.m[i] - rhs.m[i];
        }
        return result;
    }

    MAT_TYPE operator/(const MAT_TYPE& rhs) const noexcept
    {
        auto&    self = type();
        MAT_TYPE result{};
        for(size_type i = 0; i < size; ++i)
        {
            assert(rhs.m[i] > T(kEpsilon));
            result.m[i] = self.m[i] / rhs.m[i];
        }
        return result;
    }

    MAT_TYPE operator+(const_reference scalar) const noexcept
    {
        auto&    self = type();
        MAT_TYPE result{};
        for(size_type i = 0; i < size; ++i)
        {
            result.m[i] = self.m[i] + scalar;
        }
        return result;
    }

    MAT_TYPE operator-(const_reference scalar) const noexcept
    {
        auto&    self = type();
        MAT_TYPE result{};
        for(size_type i = 0; i < size; ++i)
        {
            result.m[i] = self.m[i] - scalar;
        }
        return result;
    }

    MAT_TYPE operator*(const_reference scalar) const noexcept
    {
        auto&    self = type();
        MAT_TYPE result{};
        for(size_type i = 0; i < size; ++i)
        {
            result.m[i] = self.m[i] * scalar;
        }
        return result;
    }

    MAT_TYPE operator/(const_reference scalar) const noexcept
    {
        auto&    self = type();
        MAT_TYPE result{};
        for(size_type i = 0; i < size; ++i)
        {
            assert(scalar > T(kEpsilon));
            result.m[i] = self.m[i] / scalar;
        }
        return result;
    }

    MAT_TYPE lerp(const MAT_TYPE& start, const MAT_TYPE& end, const_reference t) const noexcept
    {
        auto&    self = type();
        MAT_TYPE result{};
        for(size_type i = 0; i < size; ++i)
        {
            result.m[i] = start.m[i] + (end.m[i] - start.m[i]) * t;
        }
        return result;
    }

    friend std::ostream& operator<<(std::ostream& stream, const MAT_TYPE& matrix) noexcept
    {
        for(size_type col = 0; col < COLUMNS; ++col)
        {
            stream << "[";
            for(size_type row = 0; row < ROWS; ++row)
            {
                stream << matrix.m[col * ROWS + row];
                if(row < ROWS - 1)
                {
                    stream << "\t";
                }
            }
            stream << "]" << std::endl;
        }
        return stream;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// MAT4 SPECIFIC METHODS ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 4 && C == 4, MAT_TYPE> operator*(const MAT_TYPE& rhs) const noexcept
    {
#define M4_DOT(Ra, Cb)                                                                            \
    self.m[0 * ROWS + Ra] * rhs.m[Cb * ROWS + 0] + self.m[1 * ROWS + Ra] * rhs.m[Cb * ROWS + 1] + \
        self.m[2 * ROWS + Ra] * rhs.m[Cb * ROWS + 2] + self.m[3 * ROWS + Ra] * rhs.m[Cb * ROWS + 3];

        const auto& self = type();
        MAT_TYPE    result{};
        result.m[0]  = M4_DOT(0, 0); // col#1
        result.m[1]  = M4_DOT(1, 0);
        result.m[2]  = M4_DOT(2, 0);
        result.m[3]  = M4_DOT(3, 0);
        result.m[4]  = M4_DOT(0, 1); // col#2
        result.m[5]  = M4_DOT(1, 1);
        result.m[6]  = M4_DOT(2, 1);
        result.m[7]  = M4_DOT(3, 1);
        result.m[8]  = M4_DOT(0, 2); // col#3
        result.m[9]  = M4_DOT(1, 2);
        result.m[10] = M4_DOT(2, 2);
        result.m[11] = M4_DOT(3, 2);
        result.m[12] = M4_DOT(0, 3); // col#4
        result.m[13] = M4_DOT(1, 3);
        result.m[14] = M4_DOT(2, 3);
        result.m[15] = M4_DOT(3, 3);
        return result;
    }

    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 4 && C == 4, MAT_TYPE&> transpose() noexcept
    {
        /*
         * Consider the below layout of the matrix for calculating the cofactor matrix.
         *
         * Column-major layout for mat4
         *
         * _______________________________
         * | m[0] | m[4] | m[8]  | m[12] |
         * -------------------------------
         * | m[1] | m[5] | m[9]  | m[13] |
         * -------------------------------
         * | m[2] | m[6] | m[10] | m[14] |
         * -------------------------------
         * | m[3] | m[7] | m[11] | m[15] |
         * -------------------------------
         *
         */
        auto& self = type();
        std::swap(self.m[1], self.m[4]);
        std::swap(self.m[2], self.m[8]);
        std::swap(self.m[3], self.m[12]);
        std::swap(self.m[6], self.m[9]);
        std::swap(self.m[7], self.m[13]);
        std::swap(self.m[11], self.m[14]);
        return self;
    }

    template <size_type R = ROWS, size_type C = COLUMNS>
    [[nodiscard]] typename std::enable_if_t<R == 4 && C == 4, MAT_TYPE> transposed() const noexcept
    {
        MAT_TYPE result(type());
        return result.transpose();
    }

    template <size_type R = ROWS, size_type C = COLUMNS>
    [[nodiscard]] typename std::enable_if_t<R == 4 && C == 4, real> determinant() const noexcept
    {
        /*
         * Consider the below layout of the matrix for calculating the cofactor matrix.
         *
         * Column-major layout for mat4
         *
         * _______________________________
         * | m[0] | m[4] | m[8]  | m[12] |
         * -------------------------------
         * | m[1] | m[5] | m[9]  | m[13] |
         * -------------------------------
         * | m[2] | m[6] | m[10] | m[14] |
         * -------------------------------
         * | m[3] | m[7] | m[11] | m[15] |
         * -------------------------------
         *
         * Here we are using Laplace Expansion to calculate the determinant.
         * Determinant is given by the sum of multiplication of each element in the first column if the matrix with its
         * cofactor.
         *
         */
        value_type  result{T{0}};
        const auto& m = type().m;
        result += m[0] * (m[5] * (m[10] * m[15] - m[11] * m[14]) - m[9] * (m[6] * m[15] - m[7] * m[14]) +
                          m[13] * (m[6] * m[11] - m[7] * m[10]));
        result -= m[1] * (m[4] * (m[10] * m[15] - m[11] * m[14]) - m[8] * (m[6] * m[15] - m[7] * m[14]) +
                          m[12] * (m[6] * m[11] - m[7] * m[10]));
        result += m[2] * (m[4] * (m[9] * m[15] - m[11] * m[13]) - m[8] * (m[5] * m[15] - m[7] * m[13]) +
                          m[12] * (m[5] * m[11] - m[9] * m[7]));
        result -= m[3] * (m[4] * (m[9] * m[14] - m[10] * m[13]) - m[8] * (m[5] * m[14] - m[6] * m[13]) +
                          m[12] * (m[5] * m[10] - m[6] * m[9]));
        return result;
    }

    /*!
     * @brief Inverts this matrix. If t is not invertible, an identity matrix is returned.
     *
     * @return A reference to the inverted matrix.
     */
    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 4 && C == 4, MAT_TYPE&> invert() noexcept
    {
        /*
         * Consider the below layout of the matrix for calculating the cofactor matrix.
         *
         * Column-major layout for mat4
         *
         * _______________________________
         * | m[0] | m[4] | m[8]  | m[12] |
         * -------------------------------
         * | m[1] | m[5] | m[9]  | m[13] |
         * -------------------------------
         * | m[2] | m[6] | m[10] | m[14] |
         * -------------------------------
         * | m[3] | m[7] | m[11] | m[15] |
         * -------------------------------
         *
         */

        auto&    self = type();
        MAT_TYPE cofactor_matrix{};
        auto&    m           = self.m;
        cofactor_matrix.m[0] = (m[5] * (m[10] * m[15] - m[11] * m[14]) - m[9] * (m[6] * m[15] - m[7] * m[14]) +
                                m[13] * (m[6] * m[11] - m[7] * m[10]));
        cofactor_matrix.m[1] = -(m[4] * (m[10] * m[15] - m[11] * m[14]) - m[8] * (m[6] * m[15] - m[7] * m[14]) +
                                 m[12] * (m[6] * m[11] - m[7] * m[10]));
        cofactor_matrix.m[2] = (m[4] * (m[9] * m[15] - m[11] * m[13]) - m[8] * (m[5] * m[15] - m[7] * m[13]) +
                                m[12] * (m[5] * m[11] - m[9] * m[7]));
        cofactor_matrix.m[3] = -(m[4] * (m[9] * m[14] - m[10] * m[13]) - m[8] * (m[5] * m[14] - m[6] * m[13]) +
                                 m[12] * (m[5] * m[10] - m[6] * m[9]));
        cofactor_matrix.m[4] = -(m[1] * (m[10] * m[15] - m[11] * m[14])) - m[9] * (m[2] * m[15] - m[3] * m[14]) +
                               m[13] * (m[2] * m[11] - m[3] * m[10]);
        cofactor_matrix.m[5]  = (m[0] * (m[10] * m[15] - m[11] * m[14]) - m[8] * (m[2] * m[15] - m[3] * m[14]) +
                                m[12] * (m[2 * m[11] - m[3] * m[10]]));
        cofactor_matrix.m[6]  = -(m[0] * (m[9] * m[15] - m[11] * m[13]) - m[8] * (m[1] * m[15] - m[3] * m[13]) +
                                 m[12] * (m[1] * m[11] - m[3] * m[9]));
        cofactor_matrix.m[7]  = (m[0] * (m[9] * m[14] - m[10] * m[13]) - m[8] * (m[1] * m[14] - m[2] * m[13]) +
                                m[12] * (m[1] * m[10] - m[2] * m[9]));
        cofactor_matrix.m[8]  = (m[1] * (m[6] * m[15] - m[7] * m[14]) - m[5] * (m[2] * m[15] - m[3] * m[14]) +
                                m[13] * (m[2] * m[7] - m[3] * m[6]));
        cofactor_matrix.m[9]  = -(m[0] * (m[6] * m[15] - m[7] * m[14]) - m[4] * (m[2] * m[15] - m[3] * m[14]) +
                                 m[12] * (m[2] * m[7] - m[3] * m[6]));
        cofactor_matrix.m[10] = (m[0] * (m[5] * m[15] - m[7] * m[13]) - m[4] * (m[1] * m[15] - m[3] * m[13]) +
                                 m[12] * (m[1] * m[7] - m[3] * m[5]));
        cofactor_matrix.m[11] = -(m[0] * (m[5] * m[14] - m[6] * m[13]) - m[4] * (m[1] * m[14] - m[2] * m[13]) +
                                  m[12] * (m[1] * m[6] - m[2] * m[5]));
        cofactor_matrix.m[12] = -(m[1] * (m[6] * m[11] - m[7] * m[10]) - m[5] * (m[2] * m[11] - m[3] * m[10]) +
                                  m[9] * (m[2] * m[7] - m[3] * m[6]));
        cofactor_matrix.m[13] = (m[0] * (m[6] * m[11] - m[7] * m[10]) - m[4] * (m[2] * m[11] - m[3] * m[10]) +
                                 m[8] * (m[2] * m[7] - m[3] * m[6]));
        cofactor_matrix.m[14] = -(m[0] * (m[5] * m[11] - m[9] * m[7]) - m[4] * (m[1] * m[10] - m[9] * m[3]) +
                                  m[8] * (m[1] * m[7] - m[5] * m[3]));
        cofactor_matrix.m[15] = (m[0] * (m[5] * m[10] - m[6] * m[9]) - m[4] * (m[1] * m[10] - m[2] * m[9]) +
                                 m[8] * (m[1] * m[6] - m[2] * m[5]));

        // Calculate the determinant
        value_type determinant = m[0] * cofactor_matrix.m[0] + m[1] * cofactor_matrix.m[1] +
                                 m[2] * cofactor_matrix.m[2] + m[3] * cofactor_matrix.m[3];
        if(determinant < kEpsilon)
        {
            return MAT_TYPE{}; // return an identity matrix (maybe have a static property / method)
        }

        // Transpose the cofactor matrix to get adjugate matrix
        std::swap(cofactor_matrix.m[1], cofactor_matrix.m[4]);
        std::swap(cofactor_matrix.m[2], cofactor_matrix.m[8]);
        std::swap(cofactor_matrix.m[3], cofactor_matrix.m[12]);
        std::swap(cofactor_matrix.m[6], cofactor_matrix.m[9]);
        std::swap(cofactor_matrix.m[7], cofactor_matrix.m[13]);
        std::swap(cofactor_matrix.m[11], cofactor_matrix.m[14]);

        // Divide the adjugate matrix by the determinant to get the inverse of this matrix
        m = cofactor_matrix / determinant;
        return self;
    }

    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 4 && C == 4, MAT_TYPE> inverse() const noexcept
    {
        MAT_TYPE copy(type());
        return copy.invert();
    }

    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 4 && C == 4, MAT_TYPE> scale(const vec3& scale) const noexcept
    {
        // This is copilot generated code
        // Todo:: learn this? Why multiplying scale and not adding?
        const auto& self = type();
        MAT_TYPE    result(self);
        result.m[0] *= scale.x;
        result.m[5] *= scale.y;
        result.m[10] *= scale.z;
        return result;
    }

    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 4 && C == 4, MAT_TYPE> rotate(const vec3& axis, const real& radians) const noexcept
    {
        // This is copilot generated code
        // Todo:: learn this? how to encode rotation into a matrix
        const auto& self = type();
        MAT_TYPE    result(self);
        const real  c = real_cos(radians);
        const real  s = real_sin(radians);
        const real  t = 1 - c;
        const real  x = axis.x;
        const real  y = axis.y;
        const real  z = axis.z;
        result.m[0]   = t * x * x + c;
        result.m[1]   = t * x * y + s * z;
        result.m[2]   = t * x * z - s * y;
        result.m[4]   = t * x * y - s * z;
        result.m[5]   = t * y * y + c;
        result.m[6]   = t * y * z + s * x;
        result.m[8]   = t * x * z + s * y;
        result.m[9]   = t * y * z - s * x;
        result.m[10]  = t * z * z + c;
        return result;
    }

    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 4 && C == 4, MAT_TYPE> translate(const vec3& translation) const noexcept
    {
        const auto& self = type();
        MAT_TYPE    result(self);
        result.m[12] += translation.x;
        result.m[13] += translation.y;
        result.m[14] += translation.z;
        return result;
    }

    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 4 && C == 4, vec3> transformVector(const vec3& vector) const noexcept
    {
        const auto& self = type();
        vec3        result(vector);
        result.x = vector.x * self.m[0] + vector.y * self.m[4] + vector.z * self.m[8];
        result.y = vector.x * self.m[1] + vector.y * self.m[5] + vector.z * self.m[9];
        result.z = vector.x * self.m[2] + vector.y * self.m[6] + vector.z * self.m[10];
        return result;
    }

    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 4 && C == 4, vec3> transformPoint(const vec3& vector) const noexcept
    {
        const auto& self = type();
        vec3        result(vector);
        result.x = vector.x * self.m[0] + vector.y * self.m[4] + vector.z * self.m[8] + self.m[12];
        result.y = vector.x * self.m[1] + vector.y * self.m[5] + vector.z * self.m[9] + self.m[13];
        result.z = vector.x * self.m[2] + vector.y * self.m[6] + vector.z * self.m[10] + self.m[14];
        return result;
    }

    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 4 && C == 4, vec3> transformPoint(const vec3& vector, real& w) const noexcept
    {
        const auto& self = type();
        vec3        result(vector);
        result.x = vector.x * self.m[0] + vector.y * self.m[4] + vector.z * self.m[8] + w * self.m[12];
        result.y = vector.x * self.m[1] + vector.y * self.m[5] + vector.z * self.m[9] + w * self.m[13];
        result.z = vector.x * self.m[2] + vector.y * self.m[6] + vector.z * self.m[10] + w * self.m[14];
        w        = vector.x * self.m[3] + vector.y * self.m[7] + vector.z * self.m[11] + w * self.m[15];
        return result;
    }

    ////////////////////////////////////// MAT4 SPECIFIC METHODS END /////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////// MAT3 SPECIFIC METHODS ///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 3 && C == 3, MAT_TYPE> operator*(const MAT_TYPE& rhs) const noexcept
    {
#define M3_DOT(Ra, Cb)                                                                            \
    self.m[0 * ROWS + Ra] * rhs.m[Cb * ROWS + 0] + self.m[1 * ROWS + Ra] * rhs.m[Cb * ROWS + 1] + \
        self.m[2 * ROWS + Ra] * rhs.m[Cb * ROWS + 2];

        const auto& self = type();
        MAT_TYPE    result{};
        result.m[0] = M3_DOT(0, 0); // col#1
        result.m[1] = M3_DOT(1, 0);
        result.m[2] = M3_DOT(2, 0);
        result.m[3] = M3_DOT(0, 1); // col#2
        result.m[4] = M3_DOT(1, 1);
        result.m[5] = M3_DOT(2, 1);
        result.m[6] = M3_DOT(0, 2); // col#3
        result.m[7] = M3_DOT(1, 2);
        result.m[8] = M3_DOT(2, 2);
        return result;
    }

    template <size_type R = ROWS, size_type C = COLUMNS>
    [[nodiscard]] typename std::enable_if_t<R == 3 && C == 3, real> determinant() const noexcept
    {
        /*
         * Consider the below layout of the matrix for calculating the cofactor matrix.
         *
         * Column-major layout for mat3
         *
         * ______________________
         * | m[0] | m[3] | m[6] |
         * ----------------------
         * | m[1] | m[4] | m[7] |
         * ----------------------
         * | m[2] | m[5] | m[8] |
         * ----------------------
         *
         */
        value_type  result{T{0}};
        const auto& m = type().m;
        result += m[0] * (m[4] * m[8] - m[5] * m[7]);
        result -= m[1] * (m[3] * m[8] - m[5] * m[6]);
        result += m[2] * (m[3] * m[7] - m[4] * m[6]);
        return result;
    }

    /*!
     * @brief Inverts this matrix. If t is not invertible, an identity matrix is returned.
     *
     * @return A reference to the inverted matrix.
     */
    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 3 && C == 3, MAT_TYPE&> invert() noexcept
    {
        /*
         * Consider the below layout of the matrix for calculating the cofactor matrix.
         *
         * Column-major layout for mat3
         *
         * ______________________
         * | m[0] | m[3] | m[6] |
         * ----------------------
         * | m[1] | m[4] | m[7] |
         * ----------------------
         * | m[2] | m[5] | m[8] |
         * ----------------------
         *
         */

        auto&    self = type();
        MAT_TYPE cofactor_matrix{};
        auto&    m           = self.m;
        cofactor_matrix.m[0] = +(m[4] * m[8] - m[5] * m[7]);
        cofactor_matrix.m[1] = -(m[3] * m[8] - m[5] * m[6]);
        cofactor_matrix.m[2] = +(m[3] * m[7] - m[4] * m[6]);
        cofactor_matrix.m[3] = -(m[1] * m[8] - m[2] * m[7]);
        cofactor_matrix.m[4] = +(m[0] * m[8] - m[2] * m[6]);
        cofactor_matrix.m[5] = -(m[0] * m[7] - m[1] * m[6]);
        cofactor_matrix.m[6] = +(m[1] * m[5] - m[2] * m[4]);
        cofactor_matrix.m[7] = -(m[0] * m[5] - m[2] * m[3]);
        cofactor_matrix.m[8] = +(m[0] * m[4] - m[1] * m[3]);

        // Calculate the determinant
        value_type determinant =
            m[0] * cofactor_matrix.m[0] + m[1] * cofactor_matrix.m[1] + m[2] * cofactor_matrix.m[2];
        if(determinant < kEpsilon)
        {
            return MAT_TYPE{}; // return an identity matrix (maybe have a static property / method)
        }

        // Transpose the cofactor matrix to get adjugate matrix
        std::swap(cofactor_matrix.m[1], cofactor_matrix.m[3]);
        std::swap(cofactor_matrix.m[2], cofactor_matrix.m[6]);
        std::swap(cofactor_matrix.m[5], cofactor_matrix.m[7]);

        // Divide the adjugate matrix by the determinant to get the inverse of this matrix
        m = cofactor_matrix / determinant;
        return self;
    }

    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 3 && C == 3, MAT_TYPE> inverse() const noexcept
    {
        MAT_TYPE copy(type());
        return copy.invert();
    }

    ////////////////////////////////////// MAT3 SPECIFIC METHODS END /////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////// MAT2 SPECIFIC METHODS ///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    template <size_type R = ROWS, size_type C = COLUMNS>
    typename std::enable_if_t<R == 2 && C == 2, MAT_TYPE> operator*(const MAT_TYPE& rhs) const noexcept
    {
#define M2_DOT(Ra, Cb) self.m[0 * ROWS + Ra] * rhs.m[Cb * ROWS + 0] + self.m[1 * ROWS + Ra] * rhs.m[Cb * ROWS + 1];

        const auto& self = type();
        MAT_TYPE    result{};
        result.m[0] = M2_DOT(0, 0); // col#1
        result.m[1] = M2_DOT(1, 0);
        result.m[2] = M2_DOT(0, 1); // col#2
        result.m[3] = M2_DOT(1, 1);
        return result;
    }

    //////////////////////////////////// MAT2 SPECIFIC METHODS END ///////////////////////////////////////
};

/*
 * Consider the below layout of the matrix for calculating the cofactor matrix.
 *
 * Column-major layout for mat4
 *
 * _______________________________
 * | m[0] | m[4] | m[8]  | m[12] |
 * -------------------------------
 * | m[1] | m[5] | m[9]  | m[13] |
 * -------------------------------
 * | m[2] | m[6] | m[10] | m[14] |
 * -------------------------------
 * | m[3] | m[7] | m[11] | m[15] |
 * -------------------------------
 *
 */
struct mat4 : public TMatrix<mat4, real, 4, 4>
{
    union
    {
        struct
        {
            real m00, m01, m02, m03; // column#1
            real m10, m11, m12, m13; // column#2
            real m20, m21, m22, m23; // column#3
            real m30, m31, m32, m33; // column#4
        };

        // basis vectors notation
        struct
        {
            vec4 right;    // column#1
            vec4 up;       // column#2
            vec4 forward;  // column#3
            vec4 position; // column#4
        };

        struct
        {
            real xx, xy, xz, xw; // column#1
            real yx, yy, yz, yw; // column#2
            real zx, zy, zz, zw; // column#3
            real tx, ty, tz, tw; // column#4
        };

        //[column, row] notation
        struct
        {
            real c0r0, c0r1, c0r2, c0r3; // column#1
            real c1r0, c1r1, c1r2, c1r3; // column#2
            real c2r0, c2r1, c2r2, c2r3; // column#3
            real c3r0, c3r1, c3r2, c3r3; // column#4
        };

        //[row, column] notation
        struct
        {
            real r0c0, r1c0, r2c0, r3c0; // column#1
            real r0c1, r1c1, r2c1, r3c1; // column#2
            real r0c2, r1c2, r2c2, r3c2; // column#3
            real r0c3, r1c3, r2c3, r3c3; // column#4
        };

        // column-major notation
        struct
        {
            vec4 col0;
            vec4 col1;
            vec4 col2;
            vec4 col3;
        };

        std::array<real, 16> m;
    };

    // mat4() noexcept
    //     : right{real(1.0), 0, 0, 0}
    //     , up{0, real(1.0), 0, 0}
    //     , forward{0, 0, real(1.0), 0}
    //     , position{0, 0, 0, real(1.0)} {};

    mat4() noexcept : xx{real(1.0)}, yy{real(1.0)}, zz{real(1.0)}, tw{real(1.0)} {};

    mat4(const real* const v) noexcept
        : xx(v[0]) // col#1
        , xy(v[1])
        , xz(v[2])
        , xw(v[3])
        , yx(v[4]) // col#2
        , yy(v[5])
        , yz(v[6])
        , yw(v[7])
        , zx(v[8]) // col#3
        , zy(v[9])
        , zz(v[10])
        , zw(v[11])
        , tx(v[12]) // col#4
        , ty(v[13])
        , tz(v[14])
        , tw(v[15])
    {
    }

    mat4(const real& c0r0,
         const real& c0r1,
         const real& c0r2,
         const real& c0r3,
         const real& c1r0,
         const real& c1r1,
         const real& c1r2,
         const real& c1r3,
         const real& c2r0,
         const real& c2r1,
         const real& c2r2,
         const real& c2r3,
         const real& c3r0,
         const real& c3r1,
         const real& c3r2,
         const real& c3r3) noexcept
        : xx(c0r0) // col#1
        , xy(c0r1)
        , xz(c0r2)
        , xw(c0r3)
        , yx(c1r0) // col#2
        , yy(c1r1)
        , yz(c1r2)
        , yw(c1r3)
        , zx(c2r0) // col#3
        , zy(c2r1)
        , zz(c2r2)
        , zw(c2r3)
        , tx(c3r0) // col#4
        , ty(c3r1)
        , tz(c3r2)
        , tw(c3r3)
    {
    }

    mat4(const vec4& _right, const vec4& _up, const vec4& _forward, const vec4& _position) noexcept
        : right(_right)
        , up(_up)
        , forward(_forward)
        , position(_position)
    {
    }

    mat4& operator=(const mat4& other) noexcept = default;
    //{
    //    xx = other.xx;
    //    xy = other.xy;
    //    xz = other.xz;
    //    xw = other.xw;
    //    yx = other.yx;
    //    yy = other.yy;
    //    yz = other.yz;
    //    yw = other.yw;
    //    zx = other.zx;
    //    zy = other.zy;
    //    zz = other.zz;
    //    zw = other.zw;
    //    tx = other.tx;
    //    ty = other.ty;
    //    tz = other.tz;
    //    tw = other.tw;
    //    return *this;
    //}

    mat4& operator=(mat4&& other) noexcept = default;
};

/*
 * Consider the below layout of the matrix for calculating the cofactor matrix.
 *
 * Column-major layout for mat3
 *
 * ______________________
 * | m[0] | m[3] | m[6] |
 * ----------------------
 * | m[1] | m[4] | m[7] |
 * ----------------------
 * | m[2] | m[5] | m[8] |
 * ----------------------
 *
 */
struct mat3 : public TMatrix<mat3, real, 4, 4>
{
    union
    {
        struct
        {
            real m00, m01, m02; // column#1
            real m10, m11, m12; // column#2
            real m20, m21, m22; // column#3
        };

        // basis vectors notation
        struct
        {
            vec3 right;   // column#1
            vec3 up;      // column#2
            vec3 forward; // column#3
        };

        struct
        {
            real xx, xy, xw; // column#1
            real yx, yy, yw; // column#2
            real tx, ty, tw; // column#3
        };

        //[column, row] notation
        struct
        {
            real c0r0, c0r1, c0r2; // column#1
            real c1r0, c1r1, c1r2; // column#2
            real c2r0, c2r1, c2r2; // column#3
        };

        //[row, column] notation
        struct
        {
            real r0c0, r1c0, r2c0; // column#1
            real r0c1, r1c1, r2c1; // column#2
            real r0c2, r1c2, r2c2; // column#3
        };

        // column-major notation
        struct
        {
            vec3 col0;
            vec3 col1;
            vec3 col2;
        };

        std::array<real, 9> m;
    };

    // mat3() noexcept
    //     : right{real(1.0), 0, 0, 0}
    //     , up{0, real(1.0), 0, 0}
    //     , forward{0, 0, real(1.0), 0}
    //     , position{0, 0, 0, real(1.0)} {};

    mat3() noexcept : xx{real(1.0)}, yy{real(1.0)}, tw{real(1.0)} {};

    mat3(const real* const v) noexcept
        : xx(v[0]) // col#1
        , xy(v[1])
        , xw(v[2])
        , yx(v[4]) // col#2
        , yy(v[5])
        , yw(v[6])
        , tx(v[8]) // col#3
        , ty(v[9])
        , tw(v[10])
    {
    }

    mat3(const real& c0r0,
         const real& c0r1,
         const real& c0r2,
         const real& c1r0,
         const real& c1r1,
         const real& c1r2,
         const real& c2r0,
         const real& c2r1,
         const real& c2r2) noexcept
        : xx(c0r0) // col#1
        , xy(c0r1)
        , xw(c0r2)
        , yx(c1r0) // col#2
        , yy(c1r1)
        , yw(c1r2)
        , tx(c2r0) // col#3
        , ty(c2r1)
        , tw(c2r2)
    {
    }

    mat3(const vec3& _right, const vec3& _up, const vec3& _forward) noexcept : right(_right), up(_up), forward(_forward)
    {
    }

    mat3& operator=(const mat3& other) noexcept = default;
    //{
    //    xx = other.xx;
    //    xy = other.xy;
    //    xz = other.xz;
    //    xw = other.xw;
    //    yx = other.yx;
    //    yy = other.yy;
    //    yz = other.yz;
    //    yw = other.yw;
    //    zx = other.zx;
    //    zy = other.zy;
    //    zz = other.zz;
    //    zw = other.zw;
    //    tx = other.tx;
    //    ty = other.ty;
    //    tz = other.tz;
    //    tw = other.tw;
    //    return *this;
    //}

    mat3& operator=(mat3&& other) noexcept = default;
};

} // namespace ramanujan::experimental

#endif
