find_package(Boost COMPONENTS unit_test_framework random timer REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

include_directories(${Gps_SOURCE_DIR}/include /home/tw395/opt/include/StOpt/cdf /home/tw395/opt/include/eigen3)
add_executable(gpsTest gpsTest.cpp)
target_link_directories(gpsTest PUBLIC /home/tw395/opt/lib)
target_link_libraries(gpsTest StOpt ${Boost_RANDOM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_TIMER_LIBRARY})
