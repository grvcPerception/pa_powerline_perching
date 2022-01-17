#pragma once
#include <math_types.h>

typedef struct XInit
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
    Scalar t0;
    Scalar t1;
    Scalar t2;
    Scalar t3;

    XInit() = default;
    XInit(const std::string &path, const Scalar z_min);
    XInit(const std::string &path, const Scalar z_min, const Scalar thrust_hover);
} XInit;

typedef union {
    XInit elements;
    Scalar data[17];
} XInitVec;