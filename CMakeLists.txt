cmake_minimum_required (VERSION 2.8.7)

project (CDTGenerator)

set (CMAKE_CXX_FLAGS "-std=c++0x")

set (CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/lib/x86_64-linux-gnu)

#find_library (CGAL /usr/local/lib)

include_directories (BEFORE include)


add_subdirectory (src)


