add_executable(computeGpsCLI computeGpsCLI.cpp)

target_include_directories(computeGpsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}" "${CMAKE_SOURCE_DIR}/include/StOptCDF")

target_link_libraries(computeGpsCLI Boost::program_options Gps rapidcsv)

add_executable(permuteTraitsCLI permuteTraitsCLI.cpp)

target_include_directories(permuteTraitsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}" "${CMAKE_SOURCE_DIR}/include/StOptCDF")

target_link_libraries(permuteTraitsCLI Boost::program_options Gps rapidcsv)
