add_executable(gpsTest gpsTest.cpp)

target_include_directories(gpsTest PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}" "${StOpt_INCLUDE_DIRS}")

target_link_libraries(gpsTest Boost::random Eigen3::Eigen rapidcsv Gps)

add_executable(lwEcdfFault lwEcdfFault.cpp)

target_include_directories(lwEcdfFault PUBLIC "${Boost_INCLUDE_DIRS}" "${StOpt_INCLUDE_DIRS}")

target_link_libraries(lwEcdfFault Boost::random Eigen3::Eigen StOpt)
