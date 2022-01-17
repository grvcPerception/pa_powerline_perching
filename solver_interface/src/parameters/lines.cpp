#include <parameters/lines.h>
#include <yaml-cpp/yaml.h>

namespace YAML {
template<>
struct convert<std::array<Scalar, NLines>> {
  static Node encode(const std::array<Scalar, NLines>& rhs) {
    Node node;
    for(auto& l:rhs)
        node.push_back(l);
    return node;
  }

  static bool decode(const Node& node, std::array<Scalar, NLines>& rhs) {
    if(!node.IsSequence() || node.size() != NLines) {
      return false;
    }

    for (int i = 0; i < NLines; ++i)
        rhs[i] = node[i].as<Scalar>();
    return true;
  }
};
}

LineParameters::LineParameters(const std::string &path)
{
    YAML::Node params = YAML::LoadFile(path);

    lc_x = params["lc_x"].as<std::array<Scalar, NLines>>();
    lc_y = params["lc_y"].as<std::array<Scalar, NLines>>();
    lc_z = params["lc_z"].as<std::array<Scalar, NLines>>();
    lv_x = params["lv_x"].as<std::array<Scalar, NLines>>();
    lv_y = params["lv_y"].as<std::array<Scalar, NLines>>();
    lv_z = params["lv_z"].as<std::array<Scalar, NLines>>();
    rad_l = params["rad_l"].as<std::array<Scalar, NLines>>();
    sgm_length = params["sgm_length"].as<std::array<Scalar, NLines>>();

    for (int i = 0; i < NLines; ++i)
    {
      const Scalar norm = (Vector<3>() << lv_x[i], lv_y[i], lv_z[i]).finished().norm();
      lv_x[i] /= norm;
      lv_y[i] /= norm;
      lv_z[i] /= norm;
    }
};