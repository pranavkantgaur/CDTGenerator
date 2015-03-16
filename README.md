# CDTGenerator [![Build Status](https://travis-ci.org/pranavkantgaur/CodedPhaseShift3DScanner.svg)](https://travis-ci.org/pranavkantgaur/CodedPhaseShift3DScanner)
Implementation of the paper "Meshing piecewise linear complex using CDT" by Hang Si et. al.

# Building
The code is under testing right now but one can get fair idea of overall steps involved in implementation of Hang Si's Constrained Delaunay tetrahedralization(CDT) algorithm.

# Extension
Main objective is to extend this implementation to support adaptive Constrained Delaunay tetrahedralization for 3D domain.

# TODO:
* Use CGAL::Linear_cell_complex as PLC.

* Input should allow for more file options than ply.

* Why g++ -lgmp -lCGAL rply.cpp cdtCode.cpp does not link whereas, g++ rply.cpp cdtCode.cpp -lgmp -lCGAL does? Also there is issue in example usage of -l option given in 'man gcc' :(

* Replace naive code for finding all possible pairs of segments.

* How to remove hardcoded seeds for CGAL random number generator?

* Add termination criteria using local feature size for each vertex.

* Add functionality for Algorithm execution visualization. 


