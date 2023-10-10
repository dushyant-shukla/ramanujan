#pragma once

#include "precision.h"

namespace ramanujan::experimental
{
struct Vector3
{
    union
    {
        struct
        {
            real x;
            real y;
            real z;
        };
        real v[3];
    };

    Vector3();
    Vector3(real v);
    Vector3(real _x, real _y, real _z);
    Vector3(real* fv);

    bool operator==(const Vector3& other);
    bool operator!=(const Vector3& other);
};
}

