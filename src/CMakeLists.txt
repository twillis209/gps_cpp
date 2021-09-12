add_library(Gps gps.cpp)

target_include_directories(Gps PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}" "${StOpt_INCLUDE_DIRS}")

target_link_libraries(Gps Boost::random Boost::math_c99 Eigen3::Eigen StOpt)
