cmake_minimum_required(VERSION 3.1.0)

project(HoughTransform)

set(CMAKE_CXX_STANDARD 11)
include_directories(transform)
add_subdirectory(transform)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
    add_test(testing tests/testing)
endif(BUILD_TESTS)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Werror -g")
add_executable(HoughTransform main.cpp)

#target_link_libraries(HoughTransform transform)
