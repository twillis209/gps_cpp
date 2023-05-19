find_package(OpenMP)

# Added to fix this linking issue: https://stackoverflow.com/a/30877725
#add_definitions("-DBOOST_SYSTEM_NO_DEPRECATED")
add_definitions("-DBOOST_ERROR_CODE_HEADER_ONLY")

add_executable(computeGpsCLI computeGpsCLI.cpp)

target_include_directories(computeGpsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${CMAKE_SOURCE_DIR}/include/PPEcdf" "${Boost_INCLUDE_DIRS}")

target_link_libraries(computeGpsCLI Boost::program_options Boost::timer Gps rapidcsv OpenMP::OpenMP_CXX PPEcdf)

add_executable(permuteTraitsCLI permuteTraitsCLI.cpp)


target_include_directories(permuteTraitsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${CMAKE_SOURCE_DIR}/include/PPEcdf" "${Boost_INCLUDE_DIRS}")

target_link_libraries(permuteTraitsCLI Boost::program_options Gps rapidcsv PPEcdf)

add_executable(fitAndEvaluateEcdfsCLI fitAndEvaluateEcdfsCLI.cpp)

target_include_directories(fitAndEvaluateEcdfsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(fitAndEvaluateEcdfsCLI Boost::program_options Gps rapidcsv OpenMP::OpenMP_CXX PPEcdf)
