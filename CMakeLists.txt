cmake_minimum_required(VERSION 3.16)
project(filter_port LANGUAGES C CXX)

set(CMAKE_BUILD_TYPE "Debug")

find_package(Eigen3 REQUIRED)

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  # Enable colours on g++ 4.9 or greater
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9 OR CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 4.9)
    add_compile_options(-fdiagnostics-color=always)
  endif()
endif()

set(src main.cpp)

add_executable(filter_test ${src})

target_link_libraries(filter_test PRIVATE Eigen3::Eigen)