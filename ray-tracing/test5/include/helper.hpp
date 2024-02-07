#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

std::string datime_str();
float center_to_line_dist_sqr(std::vector<float>, std::vector<float>,
                              std::vector<float>);