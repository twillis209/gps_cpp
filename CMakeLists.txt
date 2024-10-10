# Setting this version so we can use FetchContent and to eliminate warning that <3.5 compatibility is going away
cmake_minimum_required(VERSION 3.5)

find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project(Gps VERSION 0.4)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

enable_testing()

# As per https://github.com/boostorg/cmake
FetchContent_Declare(
  Boost
  URL https://github.com/boostorg/boost/releases/download/boost-1.84.0/boost-1.84.0.tar.xz
  URL_MD5 893b5203b862eb9bbd08553e24ff146a
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

set(BOOST_INCLUDE_LIBRARIES multi_index)
FetchContent_MakeAvailable(Boost)

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
