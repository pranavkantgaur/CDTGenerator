cmake_minimum_required (VERSION 2.8.7)

project (CDTGenerator)

set (CMAKE_CXX_FLAGS "-std=c++0x")

set (CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/lib/x86_64-linux-gnu)

include_directories (BEFORE include /media/downloaded_libs/CGAL-4.6-beta1/include)

add_subdirectory (src)

# unit testing
add_subdirectory (gtest-1.7.0)
add_subdirectory (unit-tests)


