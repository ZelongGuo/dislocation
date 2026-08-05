# This is for CMakeLists.txt
# Zelong Guo, @ Potsdam, DE

# ------------------ CMake Version ------------------
cmake_minimum_required(VERSION 3.20)

# ------------------ Project Name -------------------
project(dislocs VERSION 0.1.0)
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")

# -------------------- Complier ---------------------
set(CMAKE_C_COMPILER gcc)

# ------------------ Options ------------------------
option(BUILD_LOCAL "Build for local" OFF)

# ------------------ Sub-targets --------------------
#include(${CMAKE_SOURCE_DIR}/cmake/tests.cmake)            # C programs debugging and test
if(BUILD_LOCAL)
    include(${CMAKE_SOURCE_DIR}/cmake/dislocs.cmake)       # Python C extension module for local test
else()
    include(${CMAKE_SOURCE_DIR}/cmake/dislocs_dist.cmake)  # Python C extension module for distribution
endif()
