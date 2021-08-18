# Testing library
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.6)
FetchContent_MakeAvailable(catch)
# Adds Catch2::Catch2

# Testing library
FetchContent_Declare(
  csv-parser
  GIT_REPOSITORY https://github.com/vincentlaucsb/csv-parser
  GIT_TAG 2.1.3)
FetchContent_MakeAvailable(csv-parser)

# Tests need to be added as executables first
add_executable(testGps testGps.cpp)

target_include_directories(Gps PUBLIC "${Gps_SOURCE_DIR}/include/gps"  "/home/tw395/opt/boost_1_77_0" "/home/tw395/opt/include/StOpt/cdf"  "/home/tw395/opt/include/eigen3")

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(testGps PRIVATE Gps Catch2::Catch2 csv)
# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME testGpsTest COMMAND testGps) # Command can be a target
