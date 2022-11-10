find_package(OpenMP)

add_executable(computeGpsCLI computeGpsCLI.cpp)

target_include_directories(computeGpsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${CMAKE_SOURCE_DIR}/include/PPEcdf" "${Boost_INCLUDE_DIRS}")

target_link_libraries(computeGpsCLI Boost::program_options Gps rapidcsv OpenMP::OpenMP_CXX PPEcdf)

add_executable(permuteTraitsCLI permuteTraitsCLI.cpp)

target_include_directories(permuteTraitsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(permuteTraitsCLI Boost::program_options Gps rapidcsv)

add_executable(fitAndEvaluateEcdfsCLI fitAndEvaluateEcdfsCLI.cpp)

target_include_directories(fitAndEvaluateEcdfsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(fitAndEvaluateEcdfsCLI Boost::program_options Gps rapidcsv OpenMP::OpenMP_CXX)

add_executable(testPPEcdfCLI testPPEcdfCLI.cpp)

target_include_directories(testPPEcdfCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/PPEcdf" "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(testPPEcdfCLI Boost::program_options PPEcdf rapidcsv OpenMP::OpenMP_CXX Gps)

add_executable(testIdxSort testIdxSort.cpp)

target_include_directories(testIdxSort PUBLIC "${CMAKE_SOURCE_DIR}/include/PPEcdf")

target_link_libraries(testIdxSort PPEcdf)
