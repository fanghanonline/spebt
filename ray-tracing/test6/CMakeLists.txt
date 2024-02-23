cmake_minimum_required(VERSION 3.0.0)
project(
  intersection-length
  VERSION 0.1.0
  LANGUAGES C CXX CUDA)
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
include(FetchContent)
FetchContent_Declare(
        yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG 0.8.0
)

FetchContent_GetProperties(yaml-cpp)
if(NOT yaml-cpp_POPULATED)
        message(STATUS "Fetching yaml-cpp")
        FetchContent_Populate(yaml-cpp)
        add_subdirectory(${yaml-cpp_SOURCE_DIR} ${yaml-cpp_BINARY_DIR})
endif()

include_directories(include)
include_directories(${yaml-cpp_SOURCE_DIR}/include)
add_library(helper src/helper.cpp)
add_library(kernel_functions src/kernel_functions.cu)
link_libraries(kernel_functions helper)
set_target_properties(kernel_functions PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(kernel_functions PROPERTIES CUDA_ARCHITECTURES "52;72")
add_executable(main src/main.cpp)
target_compile_options(
  kernel_functions PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: --generate-line-info
                           --use_fast_math --relocatable-device-code=true >)

target_link_libraries(main helper)
target_link_libraries(main yaml-cpp::yaml-cpp)
target_link_libraries(main kernel_functions)
set_target_properties(main PROPERTIES CUDA_ARCHITECTURES "52;72")
target_compile_options(
  main PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: --generate-line-info --use_fast_math
               --relocatable-device-code=true >)
set(CMAKE_BINARY_DIR bin)
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
