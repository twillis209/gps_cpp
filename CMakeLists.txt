cmake_minimum_required(VERSION 3.10)

project(Gps VERSION 0.1)

add_executable(gpsTest gpsTest.cpp)

configure_file(gpsConfig.h.in gpsConfig.h)

# Configured file will be written into binary tree so we need to add that to the list of paths to search for include files
target_include_directories(gpsTest PUBLIC "${PROJECT_BINARY_DIR}")
