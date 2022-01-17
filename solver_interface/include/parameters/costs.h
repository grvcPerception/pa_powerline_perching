#pragma once
#include <math_types.h>

typedef struct Costs
{
    Scalar W_t;
    Scalar W_w;
    Scalar W_xy_N;
    Scalar W_z_N;
    Scalar W_q_N;
    Scalar W_qxy_P;
    Scalar W_qz_P;
    Scalar W_v_N;
    Scalar W_w_N;
    Scalar W_pa;
    Scalar W_lc_slack;
    Scalar W_sv_slack;

    Scalar t_min;
    Scalar t_max;

    Scalar exp_decay_pa;
    Scalar exp_decay_pa_start;
    Scalar exp_decay_lc;
    Scalar exp_decay_lc_start;
    Scalar exp_decay_sv;
    Scalar exp_decay_sv_start;

    Costs() = default;
    Costs(const std::string &path);
} Costs;

typedef union {
    Costs elements;
    Scalar data[20];
} CostsVec;
