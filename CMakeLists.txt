cmake_minimum_required(VERSION 3.10)

project(Gps VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

enable_testing()

add_subdirectory(apps)
add_subdirectory(src)
add_subdirectory(test)
