#include <parameters/costs.h>
#include <yaml-cpp/yaml.h>

Costs::Costs(const std::string &path)
{
    YAML::Node params = YAML::LoadFile(path);

    W_t = params["W_t"].as<Scalar>();
    W_w = params["W_w"].as<Scalar>();
    W_xy_N = params["W_xy_N"].as<Scalar>();
    W_z_N = params["W_z_N"].as<Scalar>();
    W_q_N = params["W_q_N"].as<Scalar>();
    W_qxy_P = params["W_qxy_P"].as<Scalar>();
    W_qz_P = params["W_qz_P"].as<Scalar>();
    W_v_N = params["W_v_N"].as<Scalar>();
    W_w_N = params["W_w_N"].as<Scalar>();
    W_pa = params["W_pa"].as<Scalar>();
    W_lc_slack = params["W_lc_slack"].as<Scalar>();
    W_sv_slack = params["W_sv_slack"].as<Scalar>();
    
    t_min = params["t_min"].as<Scalar>();
    t_max = params["t_max"].as<Scalar>();

    exp_decay_pa = params["exp_decay_pa"].as<Scalar>();
    exp_decay_pa_start = params["exp_decay_pa_start"].as<Scalar>();
    exp_decay_lc = params["exp_decay_lc"].as<Scalar>();
    exp_decay_lc_start = params["exp_decay_lc_start"].as<Scalar>();
    exp_decay_sv = params["exp_decay_sv"].as<Scalar>();
    exp_decay_sv_start = params["exp_decay_sv_start"].as<Scalar>();

};