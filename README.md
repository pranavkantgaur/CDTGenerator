# CDTGenerator [![Build Status](https://travis-ci.org/pranavkantgaur/CDTGenerator.svg?branch=master)](https://travis-ci.org/pranavkantgaur/CDTGenerator) [![Join the chat at https://gitter.im/pranavkantgaur/CDTGenerator](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/pranavkantgaur/CDTGenerator?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

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

* Why ```One_dart_per_incident_cell_range<0, 2>``` produces segfault [here](http://cgal-discuss.949826.n4.nabble.com/Linear-cell-complex-Segmentation-fault-while-accessing-vertex-of-a-2-cell-td4660939.html) but ```Dart_of_orbit_range<1>``` does not?

* Template definition of ```isInfinite``` and other such functions.  

* Adding steiner point in PLC and current mesh using _edge flipping_ and combination of _edge flipping_ and _face flipping_(Refer Si's thesis).  

* Search for candidate _reference point_ only in the vertices of tetrahedrons intersecting the missing segment in function ```computeReferencePoint``` instead of all points of PLC.  

* Replace loop in current definition of ```sew2CellsFromEdge``` with _efficient_ alternative so that it can be used for larger domains as well.

* CGAL's Delaunay triangulation deals with _cospherical_ points using [symbolic perturbation](https://hal.inria.fr/inria-00166710/file/soda.pdf). How to modify the corresponding PLC accordingly?

* Input should allow for more file options than ply.

* Why g++ -lgmp -lCGAL rply.cpp cdtCode.cpp does not link whereas, g++ rply.cpp cdtCode.cpp -lgmp -lCGAL does? Also there is issue in example usage of -l option given in 'man gcc' :(

* Replace naive code for finding all possible pairs of segments.

* How to remove hardcoded seeds for CGAL random number generator?

* Add termination criteria using local feature size for each vertex.

* Add functionality for Algorithm execution visualization. 

* Add function in CGAL fork for iterating vertices of a 2-cell using ```One_dart_per_incident_cell_range<0, 2>``` in _counterclockwise_ order. It is useful while writing polygons to data file for visualization.

* Generic ```sew2CellsFromEdge``` for both ```LCC``` and ```LCCWithDartInfo```.

* How _local degeneracy removal_ alone can guarantee _uniqueness_ of Delaunay triangulation? What if a set of 5 _non-local_ points are cospherical?
