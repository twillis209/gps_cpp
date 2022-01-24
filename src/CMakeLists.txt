add_library(Gps gps.cpp)

find_package(OpenMP)

target_include_directories(Gps PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}" "${CMAKE_SOURCE_DIR}/include/StOptCDF")

target_link_libraries(Gps Eigen3::Eigen Boost::random Boost::math_c99 StOptCDF OpenMP::OpenMP_CXX)
