cmake_minimum_required(VERSION 3.0.0)
project(linesegment VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

find_package(yaml-cpp REQUIRED)
find_package(Vulkan REQUIRED)
find_package(glfw3 REQUIRED)

add_executable(main src/main.cpp)
target_link_libraries(main yaml-cpp::yaml-cpp)
target_link_libraries(main glfw Vulkan::Vulkan)
# set(CMAKE_BINARY_DIR bin)
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
