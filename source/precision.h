#ifndef RAMANUJAN_PRECISION_H
#define RAMANUJAN_PRECISION_H

#include <float.h>

namespace ramanujan
{

/*
 * Defines a real number precision. This provides us the flexibility to compile the library in single or a double
 * floating-point precision.
 */

#define SINGLE_FLOATING_POINT_PRECISION
#define LEFT_HANDED_COORDINATE_SYSTEM

#ifdef SINGLE_FLOATING_POINT_PRECISION

using real = float;

constexpr auto kPI       = 3.14159265358979323846f;
constexpr auto kPI_2     = 1.57079632679489661923f;
constexpr auto kPI_4     = 0.785398163397448309616f;
constexpr auto kEpsilon  = 0.000001f;
constexpr auto kRadToDeg = 57.2958f;
constexpr auto kDegToRad = 0.0174533f;
constexpr auto kRealMax  = FLT_MAX;

#define real_pow(x, y) powf(x, y)
#define real_sqrt(x) sqrtf(x)
#define real_acos(x) acosf(x)
#define real_sin(x) sinf(x)
#define real_cos(x) cosf(x)
#define real_abs(x) fabsf(x)
#define real_min(x, y) fminf(x, y)
#define real_max(x, y) fmaxf(x, y)
#define real_fmod(x, y) fmodf(x, y)

#else

using real = double;

constexpr auto kPI       = 3.14159265358979323846;
constexpr auto kPI_2     = 1.57079632679489661923;
constexpr auto kPI_4     = 0.785398163397448309616;
constexpr auto kEpsilon  = 0.000001;
constexpr auto kRadToDeg = 57.2958;
constexpr auto kDegToRad = 0.0174533;
constexpr auto kRealMax  = DBL_MAX;

#define real_pow(x, y) pow(x, y)
#define real_sqrt(x) sqrt(x)
#define real_acos(x) acos(x)
#define real_sin(x) sin(x)
#define real_cos(x) cos(x)
#define real_abs(x) fabs(x)
#define real_min(x, y) fmin(x, y)
#define real_max(x, y) fmax(x, y)
#define real_fmod(x, y) fmod(x, y)

#endif

} // namespace ramanujan

#endif // ! RAMANUJAN_PRECISION_H
