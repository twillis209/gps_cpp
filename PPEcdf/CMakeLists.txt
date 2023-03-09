add_library(PPEcdf PPEcdf.cpp)

find_package(OpenMP)

target_include_directories(PPEcdf PUBLIC "${CMAKE_SOURCE_DIR}/include/PPEcdf" "${Boost_INCLUDE_DIRS}" "${orbtree_SOURCE_DIR}")
