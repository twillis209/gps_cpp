add_library(Gps gps.cpp)

FetchContent_Declare(
  rapidcsv
  GIT_REPOSITORY https://github.com/d99kris/rapidcsv.git
  GIT_TAG v8.53
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(rapidcsv)

target_include_directories(Gps PUBLIC "${Gps_SOURCE_DIR}/include/gps" "${Boost_INCLUDE_DIRS}" "/opt/include/eigen-3-4.0" "/usr/local/include/StOpt")

target_link_libraries(Gps rapidcsv Boost StOpt)
