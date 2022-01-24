add_executable(computeGpsCLI computeGpsCLI.cpp)

target_include_directories(computeGpsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(computeGpsCLI Boost::program_options Gps rapidcsv)

add_executable(permuteTraitsCLI permuteTraitsCLI.cpp)

target_include_directories(permuteTraitsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(permuteTraitsCLI Boost::program_options Gps rapidcsv)

add_executable(writeFreqMapCLI writeFreqMapCLI.cpp)

target_include_directories(writeFreqMapCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(writeFreqMapCLI Boost::program_options Boost::random Boost::math_c99 Gps rapidcsv)

add_executable(computeGpsForUkbbTraitPairsCLI computeGpsForUkbbTraitPairsCLI.cpp)

target_include_directories(computeGpsForUkbbTraitPairsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(computeGpsForUkbbTraitPairsCLI Boost::program_options Gps rapidcsv)

add_executable(computeGpsForPidUkbbTraitPairsCLI computeGpsForPidUkbbTraitPairsCLI.cpp)

target_include_directories(computeGpsForPidUkbbTraitPairsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(computeGpsForPidUkbbTraitPairsCLI Boost::program_options Gps rapidcsv)
