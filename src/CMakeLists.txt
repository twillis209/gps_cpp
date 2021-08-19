add_library(Gps gps.cpp)

FetchContent_Declare(
  rapidcsv
  GIT_REPOSITORY https://github.com/d99kris/rapidcsv.git
  GIT_TAG v8.53
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(rapidcsv)

target_compile_features(Gps PRIVATE cxx_std_17)

target_include_directories(Gps PUBLIC "${Gps_SOURCE_DIR}/include/gps"  "/home/tw395/opt/boost_1_77_0" "/home/tw395/opt/include/StOpt/cdf"  "/home/tw395/opt/include/eigen3")

target_link_directories(Gps PUBLIC "/home/tw395/opt/lib" "/home/tw395/opt/boost_1_77_0/stage/lib")

target_link_libraries(Gps StOpt rapidcsv)
