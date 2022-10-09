# Testing library
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.6
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(catch)
# Adds Catch2::Catch2

# 1.3.3 builds fine
# 2.0.0-beta builds fine
# 2.0.1 does not build due to string_view-related issue
# Not sure I'm managing C++17 libraries correctly
# Testing library
#FetchContent_Declare(
#  csv-parser
#  GIT_REPOSITORY https://github.com/vincentlaucsb/csv-parser
#  GIT_TAG 2.0.0-beta
#  GIT_SHALLOW TRUE)
#FetchContent_MakeAvailable(csv-parser)

# Tests need to be added as executables first
add_executable(testGps testGps.cpp)

target_link_libraries(testGps PRIVATE Gps Catch2::Catch2 rapidcsv)

# Run with 'make test'; needs access to test data in gps/test/data
add_test(
  NAME test
  COMMAND testGps
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  )

add_executable(testFastCDF StOptCDF/testFastCDF.cpp)

#target_link_libraries(testFastCDF StOptCDF ${Boost_TIMER_LIBRARY} ${Boost_MPI_LIBRARY} ${Boost_SERIALIZATION_LIBRARY}    ${Boost_RANDOM_LIBRARY})
target_link_libraries(testFastCDF StOptCDF ${Boost_TIMER_LIBRARY} ${Boost_SERIALIZATION_LIBRARY}  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}  ${Boost_RANDOM_LIBRARY})

add_executable(testFastCDFOnSample StOptCDF/testFastCDFOnSample.cpp)

target_include_directories(testFastCDFOnSample PUBLIC "${CMAKE_SOURCE_DIR}/include/StOptCDF")

target_link_libraries(testFastCDFOnSample Eigen3::Eigen ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_TIMER_LIBRARY} StOptCDF)
