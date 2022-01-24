add_library(StOptCDF fastCDF.cpp fastCDFOnSample.cpp nDDominanceAlone.cpp)

target_include_directories(StOptCDF PUBLIC "${CMAKE_SOURCE_DIR}/include/StOptCDF")

# Eigen3 is header-only but this is still how the authors recommend the user proceed
target_link_libraries(StOptCDF Eigen3::Eigen)
