cmake_minimum_required (VERSION 2.8.7)

project (CDTGenerator)


set (CMAKE_CXX_FLAGS "-std=c++0x")

include_directories (BEFORE include)

add_subdirectory (src)


