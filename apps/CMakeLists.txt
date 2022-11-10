find_package(OpenMP)

add_executable(testIdxSort testIdxSort.cpp)

target_include_directories(testIdxSort PUBLIC "${CMAKE_SOURCE_DIR}/include/PPEcdf")

target_link_libraries(testIdxSort PPEcdf)
