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