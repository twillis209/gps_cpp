add_executable(gpsTest gpsTest.cpp)

target_include_directories(gpsTest PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}" "${CMAKE_SOURCE_DIR}/include/StOptCDF")

target_link_libraries(gpsTest Boost::random Eigen3::Eigen Gps rapidcsv)

add_executable(rgpsCLI rgpsCLI.cpp)

target_include_directories(rgpsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(rgpsCLI Gps Boost::program_options Boost::random)
