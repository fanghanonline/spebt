#include "kernel_functions.h"
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

struct vec3_dim {
  unsigned int v[3];
  unsigned int total;
  vec3_dim(unsigned int a, unsigned int b, unsigned int c) {
    v[0] = a;
    v[1] = b;
    v[2] = c;
    total = a * b * c;
  }
};
