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
add_executable(testPPEcdf testPPEcdf.cpp)

target_link_libraries(testPPEcdf PRIVATE PPEcdf Catch2::Catch2 rapidcsv)

# Run with 'make test'; needs access to test data in gps/test/data
add_test(
  NAME testPPEcdf
  COMMAND testPPEcdf
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
