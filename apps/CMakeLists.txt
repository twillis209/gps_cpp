add_executable(gpsTest gpsTest.cpp)

target_include_directories(gpsTest PUBLIC "${Gps_SOURCE_DIR}/include/gps"  "/home/tw395/opt/boost_1_77_0" "/home/tw395/opt/include/StOpt/cdf"  "/home/tw395/opt/include/eigen3")

target_link_directories(gpsTest PUBLIC "/home/tw395/opt/lib" "/home/tw395/opt/boost_1_77_0/stage/lib")

target_link_libraries(gpsTest PRIVATE Gps)
