#include <parameters/quad.h>
#include <yaml-cpp/yaml.h>
#include <iostream>

QuadParameters::QuadParameters(const std::string &path)
{
    YAML::Node params = YAML::LoadFile(path);

    g = params["g"].as<Scalar>();
    m = params["m"].as<Scalar>();
    c_t = params["c_t"].as<Scalar>();
    Jxx = params["Jxx"].as<Scalar>();
    Jyy = params["Jyy"].as<Scalar>();
    Jzz = params["Jzz"].as<Scalar>();
    lx_0 = params["lx_0"].as<Scalar>();
    lx_1 = params["lx_1"].as<Scalar>();
    lx_2 = params["lx_2"].as<Scalar>();
    lx_3 = params["lx_3"].as<Scalar>();
    ly_0 = params["ly_0"].as<Scalar>();
    ly_1 = params["ly_1"].as<Scalar>();
    ly_2 = params["ly_2"].as<Scalar>();
    ly_3 = params["ly_3"].as<Scalar>();
    max_u = params["max_u"].as<Scalar>();
    rad_x = params["rad_x"].as<Scalar>();
    rad_y = params["rad_y"].as<Scalar>();
    rad_z = params["rad_z"].as<Scalar>();
    max_rad = std::max(std::max(rad_x, rad_y), rad_z);

    const Scalar min_rad = std::min(std::min(rad_x, rad_y), rad_z);
    // For a numeric reason, this value has to be slightly different than the maximum eigenvalue of Delta
    // Simply reduce the radius by 1cm to obtain this behavior
    max_eig = std::pow(1.0 / (min_rad + 0.03 - 0.01), 2); // TODO: Set independently for each line radius

    fx = params["fx"].as<Scalar>();
    fy = params["fy"].as<Scalar>();
    pBC_x = params["pBC_x"].as<Scalar>();
    pBC_y = params["pBC_y"].as<Scalar>();
    pBC_z = params["pBC_z"].as<Scalar>();
    qBC_w = params["qBC_w"].as<Scalar>();
    qBC_x = params["qBC_x"].as<Scalar>();
    qBC_y = params["qBC_y"].as<Scalar>();
    qBC_z = params["qBC_z"].as<Scalar>();
}