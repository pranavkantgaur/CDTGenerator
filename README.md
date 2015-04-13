# CDTGenerator [![Build Status](https://travis-ci.org/pranavkantgaur/CDTGenerator.svg?branch=master)](https://travis-ci.org/pranavkantgaur/CDTGenerator)
This project is a C++ implementation of the paper "Meshing piecewise linear complex using CDT" by Hang Si et. al for generating [Constrained Delaunay tetrahedralization](http://en.wikipedia.org/wiki/Constrained_Delaunay_triangulation) for a given domain specified using set of points and _constraint_ polygons.

# Building
Following sequence of commands work on Ubuntu 14.04 LTS(64-bit):  
```mkdir build```  
```cd build```  
```cmake ..```  
```make ```  
Currently, it only builds the segment recovery module along with unit tests.  

# Dependencies
It _atleast_ works with following versions of dependencies:
* CMake(>=2.8.12) 
* CGAL(4.6)
* gtest(1.7.0)[_optional_, please see ```GTEST_ON``` switch in CMakeLists in root directory]
* gcc-4.8

# Extension
Main objective is to extend this implementation to support adaptive Constrained Delaunay tetrahedralization for 3D domain.

# TODO:

* CGAL's Delaunay triangulation deals with _cospherical_ points using [symbolic perturbation](https://hal.inria.fr/inria-00166710/file/soda.pdf). How to modify the corresponding PLC accordingly?

* Input should allow for more file options than ply.

* Why g++ -lgmp -lCGAL rply.cpp cdtCode.cpp does not link whereas, g++ rply.cpp cdtCode.cpp -lgmp -lCGAL does? Also there is issue in example usage of -l option given in 'man gcc' :(

* Replace naive code for finding all possible pairs of segments.

* How to remove hardcoded seeds for CGAL random number generator?

* Add termination criteria using local feature size for each vertex.

* Add functionality for Algorithm execution visualization. 

* Add function in CGAL fork for iterating vertices of a 2-cell using ```One_dart_per_incident_cell_range<0, 2>``` in _counterclockwise_ order. It is useful while writing polygons to data file for visualization.

