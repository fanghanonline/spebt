#include "kernel_functions.h"
#include <cstdlib>
#include <helper.hpp>

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Need config file!\n";
    exit(EXIT_FAILURE);
  }
  std::string config_fname = argv[1];
  std::cout << "Using config file: " << config_fname << "\n";
  YAML::Node config = YAML::LoadFile(config_fname);
  //   const std::string username = config["username"].as<std::string>();
  const std::vector<std::vector<float>> geomDefinition =
      config["detector geometry"].as<std::vector<std::vector<float>>>();
  std::vector<float> pointA =
      config["Debug"]["Point A"].as<std::vector<float>>();
  std::vector<float> pointB =
      config["Debug"]["Point B"].as<std::vector<float>>();
  get_solid_angle_host(geomDefinition, pointA, pointB);
  return EXIT_SUCCESS;
}