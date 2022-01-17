#include <parameters/xinit.h>
#include <yaml-cpp/yaml.h>

XInit::XInit(const std::string &path, const Scalar z_min)
{
    YAML::Node params = YAML::LoadFile(path);

    px = params["px"].as<Scalar>();
    py = params["py"].as<Scalar>();
    pz = params["pz"].as<Scalar>()-z_min;
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
    t0 = params["t0"].as<Scalar>();
    t1 = params["t1"].as<Scalar>();
    t2 = params["t2"].as<Scalar>();
    t3 = params["t3"].as<Scalar>();

    const Scalar norm = (Vector<4>() << qw, qx, qy, qz).finished().norm();
    qw /= norm;
    qx /= norm;
    qy /= norm;
    qz /= norm;
};

XInit::XInit(const std::string &path, const Scalar z_min, const Scalar thrust_hover)
{
    YAML::Node params = YAML::LoadFile(path);

    px = params["px"].as<Scalar>();
    py = params["py"].as<Scalar>();
    pz = params["pz"].as<Scalar>()-z_min;
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
    t0 = thrust_hover;
    t1 = thrust_hover;
    t2 = thrust_hover;
    t3 = thrust_hover;

    const Scalar norm = (Vector<4>() << qw, qx, qy, qz).finished().norm();
    qw /= norm;
    qx /= norm;
    qy /= norm;
    qz /= norm;
};