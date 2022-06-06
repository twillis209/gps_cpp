find_package(OpenMP)

add_executable(computeGpsCLI computeGpsCLI.cpp)

target_include_directories(computeGpsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(computeGpsCLI Boost::program_options Gps rapidcsv OpenMP::OpenMP_CXX DataFrame)

add_executable(computeGpsNoPerturbCLI computeGpsNoPerturbCLI.cpp)

target_include_directories(computeGpsNoPerturbCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(computeGpsNoPerturbCLI Boost::program_options Gps rapidcsv)

add_executable(permuteTraitsCLI permuteTraitsCLI.cpp)

target_include_directories(permuteTraitsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(permuteTraitsCLI Boost::program_options Gps rapidcsv)

add_executable(computeGpsForUkbbTraitPairsCLI computeGpsForUkbbTraitPairsCLI.cpp)

target_include_directories(computeGpsForUkbbTraitPairsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(computeGpsForUkbbTraitPairsCLI Boost::program_options Gps rapidcsv)

add_executable(computeGpsForPidUkbbTraitPairsCLI computeGpsForPidUkbbTraitPairsCLI.cpp)

target_include_directories(computeGpsForPidUkbbTraitPairsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(computeGpsForPidUkbbTraitPairsCLI Boost::program_options Gps rapidcsv)

add_executable(fitAndEvaluateEcdfsCLI fitAndEvaluateEcdfsCLI.cpp)

target_include_directories(fitAndEvaluateEcdfsCLI PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(fitAndEvaluateEcdfsCLI Boost::program_options Gps rapidcsv)

add_executable(permuteProblem permuteProblem.cpp)

target_include_directories(permuteProblem PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(permuteProblem Boost::program_options Boost::random Gps rapidcsv)

add_executable(perturbPvalues perturbPvalues.cpp)

target_include_directories(perturbPvalues PUBLIC "${CMAKE_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}")

target_link_libraries(perturbPvalues Boost::program_options Boost::random Gps rapidcsv)
