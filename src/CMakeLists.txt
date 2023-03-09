add_library(Gps gps.cpp)

find_package(OpenMP)

target_include_directories(Gps PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(Gps Boost::random Boost::math_c99 OpenMP::OpenMP_CXX)
