#include "kernel_functions.h"
#include <cstdlib>
#include <helper.hpp>

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    std::cout << "Need config file!\n";
    exit(EXIT_FAILURE);
  }
  std::string config_fname = argv[1];
  std::cout << "Using config file: " << config_fname << "\n";
  YAML::Node config = YAML::LoadFile(config_fname);
  //   const std::string username = config["username"].as<std::string>();

  get_result_host(config);
  return EXIT_SUCCESS;
}