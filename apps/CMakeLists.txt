add_executable(gpsTest gpsTest.cpp)

target_include_directories(gpsTest PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}" "${CMAKE_SOURCE_DIR}/include/StOptCDF")

target_link_libraries(gpsTest Boost::random Eigen3::Eigen Gps)

add_executable(rgpsCLI rgpsCLI.cpp)

target_include_directories(rgpsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(rgpsCLI Gps Boost::program_options Boost::random)

#add_executable(testPermutation testPermutation.cpp)

#target_include_directories(testPermutation PUBLIC "${CMAKE_SOURCE_DIR}/include/gps")

#target_link_libraries(testPermutation Gps rapidcsv)

add_executable(lwIssue lwIssue.cpp)

target_include_directories(lwIssue PUBLIC "/home/tw395/opt/include/StOpt/cdf")

target_link_directories(lwIssue PUBLIC "/home/tw395/opt/lib")

target_link_libraries(lwIssue Eigen3::Eigen Boost::random StOpt)
