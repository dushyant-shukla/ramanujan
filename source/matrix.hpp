#ifndef RAMANUJAN_MATRIX
#define RAMANUJAN_MATRIX

#include "constants.h"
#include "precision.h"
#include "vector.hpp"

#include <array>

namespace ramanujan::experimental
{

template <typename T, std::size_t ROWS, std::size_t COLUMNS>
class TMatrix
{
public:
    using value_type      = T;
    using size_type       = std::size_t;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using pointer         = value_type*;
    using const_pointer   = const value_type*;

    constexpr TMatrix() noexcept                          = default;
    constexpr TMatrix(const TMatrix&) noexcept            = default;
    constexpr TMatrix(TMatrix&&) noexcept                 = default;
    constexpr TMatrix& operator=(const TMatrix&) noexcept = default;
    constexpr TMatrix& operator=(TMatrix&&) noexcept      = default;

    constexpr reference operator()(size_type row, size_type column) noexcept { return m_data[column * ROWS + row]; }

    constexpr const_reference operator()(size_type row, size_type column) const noexcept
    {
        return m_data[column * ROWS + row];
    }

    // constexpr reference operator[](size_type index) noexcept { return data_[index]; }

    // constexpr const_reference operator[](size_type index) const noexcept { return data_[index]; }

    constexpr size_type rows() const noexcept { return ROWS; }

    constexpr size_type columns() const noexcept { return COLUMNS; }

    constexpr size_type size() const noexcept { return ROWS * COLUMNS; }

    constexpr pointer data() const noexcept { return data_.data(); }

    constexpr const_pointer data() noexcept { return data_.data(); }

protected:
    std::array<T, ROWS * COLUMNS> m_data;
};

template <typename T>
struct TMatrix3 : public TMatrix<T, 3, 3>
{
};

template <typename T>
struct TMatrix4 : public TMatrix<T, 4, 4>
{
};

using mat3 = TMatrix3<real>;
using mat4 = TMatrix4<real>;

} // namespace ramanujan::experimental

#endif
