cmake_minimum_required (VERSION 2.8.7)

project (CDTGenerator)

IF (UNIX)
	set (CMAKE_CXX_FLAGS "-std=c++0x -g")
	set (CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)
	include_directories (BEFORE include)

ELSEIF (MSVC)
	#set (CMAKE_CXX_FLAGS = "/Zi")
	set (CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} )
	#set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "C:\\Program Files\\CGAL-4.6.1\\include\\")
	include_directories (BEFORE include "C:\\Program Files\\CGAL-4.6.1\\include\\" "C:\\Program Files\\boost_1_58_0\\")
ENDIF()	


add_subdirectory (src)

# unit testing, ON/OFF unit tests
set (GTEST_ON false)

if (GTEST_ON)
	add_subdirectory (unit-tests)
endif()

