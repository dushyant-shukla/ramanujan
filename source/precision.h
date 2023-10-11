#pragma once

namespace ramanujan
{

#define SINGLE_FLOATING_POINT_PRECISION

#ifdef SINGLE_FLOATING_POINT_PRECISION

/*
 * Defines a real number precision. This provides us the flexibility to compile the library in single or a double
 * floating-point precision.
 */
using real = float;

#define real_pow(x, y) powf(x, y)

#else

using real = double;
#endif

} // namespace ramanujan
