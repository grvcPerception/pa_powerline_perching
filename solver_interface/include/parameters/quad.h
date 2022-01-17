#pragma once
#include "math_types.h"

typedef struct QuadParameters
{
    Scalar g;
    Scalar m;
    Scalar c_t;
    Scalar Jxx;
    Scalar Jyy;
    Scalar Jzz;
    Scalar lx_0;
    Scalar ly_0;
    Scalar lx_1;
    Scalar ly_1;
    Scalar lx_2;
    Scalar ly_2;
    Scalar lx_3;
    Scalar ly_3;
    Scalar max_u;
    Scalar rad_x;
    Scalar rad_y;
    Scalar rad_z;
    Scalar max_rad;
    Scalar max_eig;
    Scalar fx;
    Scalar fy;
    Scalar pBC_x;
    Scalar pBC_y;
    Scalar pBC_z;
    Scalar qBC_w;
    Scalar qBC_x;
    Scalar qBC_y;
    Scalar qBC_z;

    QuadParameters(const std::string &path);
    QuadParameters() = delete;
} QuadParameters;

typedef union {
    QuadParameters elements;
    Scalar data[29];
} QuadParametersVec;