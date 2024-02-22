find_package(OpenMP)

# Added to fix this linking issue: https://stackoverflow.com/a/30877725
#add_definitions("-DBOOST_SYSTEM_NO_DEPRECATED")
add_definitions("-DBOOST_ERROR_CODE_HEADER_ONLY")

add_executable(computeGpsCLI computeGpsCLI.cpp)

target_include_directories(computeGpsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${CMAKE_SOURCE_DIR}/include/PPEcdf" "${Boost_INCLUDE_DIRS}")

target_link_libraries(computeGpsCLI PUBLIC Gps rapidcsv OpenMP::OpenMP_CXX PPEcdf )
target_link_libraries(computeGpsCLI PRIVATE CLI11::CLI11)

add_executable(permuteTraitsCLI permuteTraitsCLI.cpp)

target_include_directories(permuteTraitsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${CMAKE_SOURCE_DIR}/include/PPEcdf" "${Boost_INCLUDE_DIRS}")

target_link_libraries(permuteTraitsCLI PUBLIC Gps rapidcsv PPEcdf)
target_link_libraries(permuteTraitsCLI PRIVATE CLI11::CLI11)

add_executable(fitAndEvaluateEcdfsCLI fitAndEvaluateEcdfsCLI.cpp)

target_include_directories(fitAndEvaluateEcdfsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(fitAndEvaluateEcdfsCLI PUBLIC Gps rapidcsv OpenMP::OpenMP_CXX PPEcdf)
target_link_libraries(fitAndEvaluateEcdfsCLI PRIVATE CLI11::CLI11)
