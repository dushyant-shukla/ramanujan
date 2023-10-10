#pragma once

namespace ramanujan
{

#if 1

/*
 * Defines a real number precision. This provides us the flexibility to compile the library in single or a double
 * floating-point precision.
 */
using real = float;

real pow(real x, real y) {}

#define real_pow(x, y) powf(x, y)

#else

#endif

} // namespace ramanujan
