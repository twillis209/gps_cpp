# Testing library
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.6
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(catch)
# Adds Catch2::Catch2

# Tests need to be added as executables first
add_executable(testGps testGps.cpp)

target_link_libraries(testGps PRIVATE PPEcdf Gps Catch2::Catch2 rapidcsv)

# Run with 'make test'; needs access to test data in gps/test/data
add_test(
  NAME test
  COMMAND testGps
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  )
