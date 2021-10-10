add_executable(gpsTest gpsTest.cpp)

target_include_directories(gpsTest PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}" "${CMAKE_SOURCE_DIR}/include/StOptCDF")

target_link_libraries(gpsTest Boost::random Eigen3::Eigen Gps rapidcsv)

add_executable(rgpsCLI rgpsCLI.cpp)

target_include_directories(rgpsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(rgpsCLI Gps Boost::program_options Boost::random)

add_executable(permuteGpsCLI permuteGpsCLI.cpp)

target_include_directories(permuteGpsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(permuteGpsCLI Gps Boost::program_options Boost::random rapidcsv)

add_executable(printCdfs printCdfs.cpp)

target_include_directories(printCdfs PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(printCdfs Gps rapidcsv)

add_executable(computeGpsCLI computeGpsCLI.cpp)

target_include_directories(computeGpsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}" "${CMAKE_SOURCE_DIR}/include/StOptCDF")

target_link_libraries(computeGpsCLI Boost::random Boost::program_options Gps rapidcsv)

add_executable(permuteUkbbTraitsCLI permuteUkbbTraitsCLI.cpp)

target_include_directories(permuteUkbbTraitsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}" "${CMAKE_SOURCE_DIR}/include/StOptCDF")

target_link_libraries(permuteUkbbTraitsCLI Boost::random Boost::program_options Gps rapidcsv)
