cmake_minimum_required(VERSION 3.1.0)

project(HoughTransform)

set(CMAKE_CXX_STANDARD 11)
include_directories(transform)
add_subdirectory(transform)

if(BUILD_TESTS)
    add_subdirectory(tests)
endif(BUILD_TESTS)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Werror")
add_executable(HoughTransform main.cpp)

#target_link_libraries(HoughTransform transform)
