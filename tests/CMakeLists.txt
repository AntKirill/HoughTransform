cmake_minimum_required(VERSION 3.1.0)
project(tests)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-pthread")

include_directories(../transform)
include_directories(gtest)
add_subdirectory(gtest)

add_executable(testing "testing.cpp")
target_link_libraries(testing googleTests)
