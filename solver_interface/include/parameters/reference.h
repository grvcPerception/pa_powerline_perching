#pragma once
#include <math_types.h>

typedef struct Reference
{
    Scalar px;
    Scalar py;
    Scalar pz;
    Scalar qw;
    Scalar qx;
    Scalar qy;
    Scalar qz;
    Scalar vx;
    Scalar vy;
    Scalar vz;
    Scalar wx;
    Scalar wy;
    Scalar wz;

    Reference() = default;
    Reference(const std::string &path);
} Reference;

typedef union {
    Reference elements;
    Scalar data[13];
} ReferenceVec;