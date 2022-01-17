#include <parameters/reference.h>
#include <yaml-cpp/yaml.h>

Reference::Reference(const std::string &path)
{
    YAML::Node params = YAML::LoadFile(path);

    px = params["px"].as<Scalar>();
    py = params["py"].as<Scalar>();
    pz = params["pz"].as<Scalar>();
    qw = params["qw"].as<Scalar>();
    qx = params["qx"].as<Scalar>();
    qy = params["qy"].as<Scalar>();
    qz = params["qz"].as<Scalar>();
    vx = params["vx"].as<Scalar>();
    vy = params["vy"].as<Scalar>();
    vz = params["vz"].as<Scalar>();
    wx = params["wx"].as<Scalar>();
    wy = params["wy"].as<Scalar>();
    wz = params["wz"].as<Scalar>();

    const Scalar norm = (Vector<4>() << qw, qx, qy, qz).finished().norm();
    qw /= norm;
    qx /= norm;
    qy /= norm;
    qz /= norm;
};