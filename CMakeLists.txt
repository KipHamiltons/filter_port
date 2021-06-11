cmake_minimum_required(VERSION 3.16)
project(filter_port LANGUAGES C CXX)

set(CMAKE_BUILD_TYPE "Debug")

set(CMAKE_CXX_STANDARD 20)

option(ENABLE_CLANG_TIDY "Enable building with clang-tidy checks." OFF)
if(ENABLE_CLANG_TIDY)
  find_package(PythonInterp 3 REQUIRED)
  set(CMAKE_CXX_CLANG_TIDY "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/clang-tidy.py"
                           "${PROJECT_BINARY_DIR}/clang-tidy-fixes" clang-tidy
  )
  set(CMAKE_C_CLANG_TIDY ${CMAKE_CXX_CLANG_TIDY})

  # Create a target that will apply clang-tidy fixes to the codebase
  add_custom_target(
    apply-clang-tidy
    COMMAND clang-apply-replacements --format --style=file --style-config="${PROJECT_SOURCE_DIR}"
            "${PROJECT_BINARY_DIR}/clang-tidy-fixes"
    COMMENT "Applying fixes from clang-tidy to the codebase."
  )
endif()

find_package(Eigen3 REQUIRED)

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  # Enable colours on g++ 4.9 or greater
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9 OR CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 4.9)
    add_compile_options(-fdiagnostics-color=always)
  endif()
endif()

set(src main.cpp kalman.cpp matrix.cpp orientation.cpp approximations.cpp tasks.hpp)

add_executable(filter_test ${src})

target_link_libraries(filter_test PRIVATE Eigen3::Eigen)