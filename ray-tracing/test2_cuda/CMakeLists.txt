cmake_minimum_required(VERSION 3.0.0)
# project(
#   linesegment
#   VERSION 0.1.0
#   LANGUAGES C CXX CUDA)
project(
  linesegment
  VERSION 0.1.0
  LANGUAGES C CXX)
if(NOT DEFINED CMAKE_CUDA_STANDARD)
  set(CMAKE_CUDA_STANDARD 11)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(CTest)
enable_testing()
enable_language(CUDA)
include(CheckLanguage)
check_language(CUDA)
find_package(yaml-cpp REQUIRED)

include_directories(include)

add_library(helper src/helper.cpp)
add_library(kernel_functions src/kernel_functions.cu)
link_libraries(kernel_functions helper)
set_target_properties(kernel_functions PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(kernel_functions PROPERTIES CUDA_ARCHITECTURES "52;72")
add_executable(main src/main.cpp)
target_link_libraries(main helper)
target_link_libraries(main yaml-cpp::yaml-cpp)
target_link_libraries(main kernel_functions)
set_target_properties(main PROPERTIES CUDA_ARCHITECTURES "52;72")
set(CMAKE_BINARY_DIR bin)
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
