# g++ -std=c++0x -g src/rply.cpp src/cdtGen.cpp -I -lgmp -lmpfr -lCGAL -o output

cmake_minimum_required (VERSION 2.8.11)

project (CDTGenerator)


set (CMAKE_CXX_FLAGS "-std=c++0x")

include_directories (BEFORE include)

add_subdirectory (src)


