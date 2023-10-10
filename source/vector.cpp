#include "vector.h"
#include "constants.h"

namespace ramanujan::experimental
{

///////////////////////// Vector3 begins /////////////////////////////////

Vector3::Vector3(): x(0.0f), y(0.0f), z(0.0f)
{
    v[0] = 0.0f;
    v[1] = 0.0f;
    v[2] = 0.0f;
}

Vector3::Vector3(real _v) : x(_v), y(_v), z(_v)
{
    v[0] = _v;
    v[1] = _v;
    v[2] = _v;
}

Vector3::Vector3(real _x, real _y, real _z) : x(_x), y(_y), z(_z)
{
    v[0] = _x;
    v[1] = _y;
    v[2] = _z;
}

Vector3::Vector3(real* fv) : x(fv[0]), y(fv[1]), z(fv[2])
{
    v[0] = fv[0];
    v[1] = fv[1];
    v[2] = fv[2];
}

bool Vector3::operator==(const Vector3& other)
{
    Vector3 diff(*this - other);
    return LengthSq(diff) < ramanujan::Constants::EPSILON;
}

bool Vector3::operator!=(const Vector3& other)
{
    return !(*this == other);
}

///////////////////////// Vector3 ends /////////////////////////////////

}