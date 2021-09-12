add_executable(gpsTest gpsTest.cpp)

target_include_directories(gpsTest PUBLIC "${Gps_SOURCE_DIR}/include/gps" ${Boost_INCLUDE_DIRS} "/opt/include/eigen-3-4.0/Eigen" "/usr/local/include/StOpt")

target_link_libraries(gpsTest PRIVATE Gps)
