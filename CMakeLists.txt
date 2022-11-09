cmake_minimum_required(VERSION 3.10)

# TODO remove if exporting
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project(Gps VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

enable_testing()

# NB: Only need to list compiled libraries here
find_package(Boost 1.77 COMPONENTS random math_c99 program_options timer serialization unit_test_framework REQUIRED)
find_package(Eigen3 3.3 REQUIRED NAMES Eigen3 Eigen eigen NO_MODULE)

FetchContent_Declare(
  rapidcsv
  GIT_REPOSITORY https://github.com/d99kris/rapidcsv.git
  GIT_TAG v8.53
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(rapidcsv)

add_subdirectory(StOptCDF)
add_subdirectory(src)
add_subdirectory(apps)
add_subdirectory(test)
add_subdirectory(PPEcdf)
