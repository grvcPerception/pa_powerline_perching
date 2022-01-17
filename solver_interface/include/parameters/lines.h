#pragma once
#include <math_types.h>

#define NLines 3

typedef struct LineParameters
{
    std::array<Scalar,NLines> lc_x;
    std::array<Scalar,NLines> lc_y;
    std::array<Scalar,NLines> lc_z;
    std::array<Scalar,NLines> lv_x;
    std::array<Scalar,NLines> lv_y;
    std::array<Scalar,NLines> lv_z;
    std::array<Scalar,NLines> rad_l;
    std::array<Scalar,NLines> sgm_length;

    LineParameters(const std::string &path);
    LineParameters() = delete;
} LineParameters;

typedef union {
    LineParameters elements;
    Scalar data[NLines*8];
} LineParametersVec;