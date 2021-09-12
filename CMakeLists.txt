cmake_minimum_required(VERSION 3.10)

# TODO remove if exporting
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project(Gps VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

enable_testing()

find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED NAMES Eigen eigen eigen-3.4.0)
find_package(StOpt REQUIRED NAMES STOPT StOpt)

FetchContent_Declare(
  rapidcsv
  GIT_REPOSITORY https://github.com/d99kris/rapidcsv.git
  GIT_TAG v8.53
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(rapidcsv)

add_subdirectory(apps)
add_subdirectory(src)
add_subdirectory(StOptCDF)
add_subdirectory(test)
