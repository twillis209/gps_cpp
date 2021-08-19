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
FetchContent_Declare(
  csv-parser
  GIT_REPOSITORY https://github.com/vincentlaucsb/csv-parser
  GIT_TAG 2.0.0-beta
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(csv-parser)

# Tests need to be added as executables first
add_executable(testGps testGps.cpp)

target_include_directories(testGps PRIVATE "${csv-parser_SOURCE_DIR}/single_include")

# csv-parser recommends this link, but it's header-only, no?
#target_link_libraries(testGps PRIVATE Gps Catch2::Catch2 csv)

target_link_libraries(testGps PRIVATE Gps Catch2::Catch2)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(
  NAME test
  COMMAND testGps
  WORKING_DIRECTORY ${Gps_ROOT_DIR}
)
