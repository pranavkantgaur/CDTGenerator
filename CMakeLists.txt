cmake_minimum_required (VERSION 2.8.7)

project (CDTGenerator)

set (CMAKE_CXX_FLAGS "-std=c++0x")

set (CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)

include_directories (BEFORE include)

add_subdirectory (src)

# unit testing, ON/OFF unit tests
set (GTEST_ON true)

if (GTEST_ON)
	add_subdirectory (unit-tests)
endif()

