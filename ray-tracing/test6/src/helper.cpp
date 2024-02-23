#include "helper.hpp"
std::string datime_str() {
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);

  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y%m%d%H%M%S");
  auto str = oss.str();

  //   std::cout << str << std::endl;
  return std::string(str);
}

float center_to_line_dist_sqr(std::vector<float> cuboid,
                              std::vector<float> pointA,
                              std::vector<float> pointB) {

  // Vector OA, OB
  // O is the center of the cuboid
  float oa[3], ob[3];
  float ab_sqr = 0;
  for (int i : {0, 1, 2}) {
    float component_i = 0.5 * (cuboid[2 * i] + cuboid[2 * i + 1]);
    oa[i] = component_i - pointA[i];
    ob[i] = component_i - pointB[i];
    ab_sqr += (pointB[i] - pointA[i]) * (pointB[i] - pointA[i]);
  }
  // OA X OB coss product x, y, z components
  float cpx = oa[1] * ob[2] - oa[2] * ob[1];
  float cpy = oa[2] * ob[0] - oa[0] * ob[2];
  float cpz = oa[0] * ob[1] - oa[1] * ob[0];
  // std::cout << "d_AB sqr=" << ab_sqr << "\n";
  return (cpx * cpx + cpy * cpy + cpz * cpz) / ab_sqr;
}