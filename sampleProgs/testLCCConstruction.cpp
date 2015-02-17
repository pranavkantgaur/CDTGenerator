/* 
 * Minimal code using CGAL's LCC package for constructing linear cell complex of a cube 
 */
#include <iostream>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Linear_cell_complex.h>

using namespace std;
using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;
//typedef Linear_cell_complex_traits<3> LCCTraits;
typedef Linear_cell_complex<3, 3> LCC;
typedef Point_3<K> Point;

int main()
{
	vector<Point> cubePoints;

	// create storage for points of a cube
	cubePoints.push_back(Point(0, 0, 0)); // points of a cube
	cubePoints.push_back(Point(0, 0, 1));
	cubePoints.push_back(Point(0, 1, 0));
	cubePoints.push_back(Point(0, 1, 1));
	cubePoints.push_back(Point(1, 0, 0));
	cubePoints.push_back(Point(1, 0, 1));
	cubePoints.push_back(Point(1, 1, 0));
	cubePoints.push_back(Point(1, 1, 1));

	LCC cubeLCC;

	// add triangles representing cube to lcc
	cubeLCC.make_triangle(cubePoints[0], cubePoints[4], cubePoints[6]);
	cubeLCC.make_triangle(cubePoints[6], cubePoints[2], cubePoints[0]);
	cubeLCC.make_triangle(cubePoints[4], cubePoints[5], cubePoints[7]);
	cubeLCC.make_triangle(cubePoints[7], cubePoints[6], cubePoints[4]);
	cubeLCC.make_triangle(cubePoints[5], cubePoints[1], cubePoints[3]);
	cubeLCC.make_triangle(cubePoints[3], cubePoints[7], cubePoints[5]);
	cubeLCC.make_triangle(cubePoints[2], cubePoints[3], cubePoints[1]);
	cubeLCC.make_triangle(cubePoints[1], cubePoints[0], cubePoints[2]);
	cubeLCC.make_triangle(cubePoints[6], cubePoints[7], cubePoints[3]);
	cubeLCC.make_triangle(cubePoints[3], cubePoints[2], cubePoints[6]);
	cubeLCC.make_triangle(cubePoints[0], cubePoints[1], cubePoints[5]);
	cubeLCC.make_triangle(cubePoints[5], cubePoints[4], cubePoints[0]);


	// Verify correctness
	unsigned int nVertices = 0;
	for (LCC::One_dart_per_cell_range<0>::iterator vIter = cubeLCC.one_dart_per_cell<0>().begin(); vIter != cubeLCC.one_dart_per_cell<0>().end(); vIter++)
		nVertices++;

	cout << "Number of vertices: " << nVertices << "\n";	

	
	unsigned int nSegments = 0;
	for (LCC::One_dart_per_cell_range<1>::iterator sIter = cubeLCC.one_dart_per_cell<1>().begin(); sIter != cubeLCC.one_dart_per_cell<1>().end(); sIter++)
		nSegments++;

	cout << "Number of segments: " << nSegments << "\n";	


	unsigned int nFaces = 0;
	for (LCC::One_dart_per_cell_range<2>::iterator fIter = cubeLCC.one_dart_per_cell<2>().begin(); fIter != cubeLCC.one_dart_per_cell<2>().end(); fIter++)
		nFaces++;

	cout << "Number of faces: " << nFaces << "\n";	

	return 0;
}

