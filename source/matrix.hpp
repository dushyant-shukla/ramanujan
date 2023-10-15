#ifndef RAMANUJAN_VECTOR
#define RAMANUJAN_VECTOR

#include "constants.h"
#include "precision.h"

#include <array>

namespace ramanujan::experimental
{

template <typename T, std::size_t ROWS, std::size_t COLUMNS>
class Matrix
{
public:
    using value_type      = T;
    using size_type       = std::size_t;
    using reference       = value_type&;
    using const_reference = const value_type&;

    constexpr Matrix() noexcept = default;

    constexpr Matrix(const Matrix&) noexcept = default;
    constexpr Matrix(Matrix&&) noexcept      = default;

    constexpr Matrix& operator=(const Matrix&) noexcept = default;
    constexpr Matrix& operator=(Matrix&&) noexcept      = default;

    constexpr Matrix(const std::array<T, ROWS * COLUMNS>& data) noexcept : data_(data) {}

    constexpr Matrix(std::array<T, ROWS * COLUMNS>&& data) noexcept : data_(std::move(data)) {}

    constexpr reference operator()(size_type row, size_type column) noexcept { return data_[row * COLUMNS + column]; }

    constexpr const_reference operator()(size_type row, size_type column) const noexcept
    {
        return data_[row * COLUMNS + column];
    }

    constexpr reference operator[](size_type index) noexcept { return data_[index]; }

    constexpr const_reference operator[](size_type index) const noexcept { return data_[index]; }

    constexpr size_type rows() const noexcept { return ROWS; }

    constexpr size_type columns() const noexcept { return COLUMNS; }

    constexpr size_type size() const noexcept { return ROWS * COLUMNS; }

    constexpr const T* data() const noexcept { return data_.data(); }

    constexpr T* data() noexcept { return data_.data(); }
};

} // namespace ramanujan::experimental

#endif
