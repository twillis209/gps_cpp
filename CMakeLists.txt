# Setting this version so we can use FetchContent
cmake_minimum_required(VERSION 3.14)

find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project(Gps VERSION 0.3)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

enable_testing()

find_package(Boost 1.59)

FetchContent_Declare(
  rapidcsv
  GIT_REPOSITORY https://github.com/d99kris/rapidcsv.git
  GIT_TAG v8.53
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(rapidcsv)

FetchContent_Declare(
  cli11_proj
  QUIET
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.3.2
)
FetchContent_MakeAvailable(cli11_proj)


add_subdirectory(src)
add_subdirectory(apps)
add_subdirectory(test)
add_subdirectory(PPEcdf)
