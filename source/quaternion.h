#pragma once

#include "vector3.h"

namespace ramanujan
{

struct Quaternion
{
    union
    {
        struct
        {
            float x;
            float y;
            float z;
            float w;
        };

        struct
        {
            Vector3 vector;
            float   scalar;
        };

        float v[4];
    };
};

} // namespace ramanujan