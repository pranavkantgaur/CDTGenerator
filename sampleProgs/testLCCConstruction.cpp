/* 
 * Minimal code using CGAL's LCC package for constructing linear cell complex of a cube 
 */
#include <iostream>
#include <vector>
#include <unordered_set>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Linear_cell_complex.h>

using namespace std;
using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;
typedef Linear_cell_complex<3, 3> LCC;
typedef Point_3<K> Point;
typedef LCC::Dart_handle DartHandle;


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
	DartHandle d1 = cubeLCC.make_triangle(cubePoints[0], cubePoints[4], cubePoints[6]);
	DartHandle d2 = cubeLCC.make_triangle(cubePoints[6], cubePoints[2], cubePoints[0]);
	DartHandle d3 = cubeLCC.make_triangle(cubePoints[4], cubePoints[5], cubePoints[7]);
	DartHandle d4 = cubeLCC.make_triangle(cubePoints[7], cubePoints[6], cubePoints[4]);
	DartHandle d5 = cubeLCC.make_triangle(cubePoints[5], cubePoints[1], cubePoints[3]);
	DartHandle d6 = cubeLCC.make_triangle(cubePoints[3], cubePoints[7], cubePoints[5]);
	DartHandle d7 = cubeLCC.make_triangle(cubePoints[2], cubePoints[3], cubePoints[1]);
	DartHandle d8 = cubeLCC.make_triangle(cubePoints[1], cubePoints[0], cubePoints[2]);
	DartHandle d9 = cubeLCC.make_triangle(cubePoints[6], cubePoints[7], cubePoints[3]);
	DartHandle d10 = cubeLCC.make_triangle(cubePoints[3], cubePoints[2], cubePoints[6]);
	DartHandle d11 = cubeLCC.make_triangle(cubePoints[0], cubePoints[1], cubePoints[5]);
	DartHandle d12 = cubeLCC.make_triangle(cubePoints[5], cubePoints[4], cubePoints[0]);


	// Gluing...
	cubeLCC.sew<2>(d1, d2);
	cubeLCC.sew<2>(d3, d4);
	cubeLCC.sew<2>(d5, d6);
	cubeLCC.sew<2>(d7, d8);
	cubeLCC.sew<2>(d9, d10);
	cubeLCC.sew<2>(d11, d12);

	cubeLCC.sew<2>(cubeLCC.beta(d1, 1), cubeLCC.beta(d4, 1));	
	cubeLCC.sew<2>(cubeLCC.beta(d1, 1, 1), cubeLCC.beta(d12, 1)); 
	cubeLCC.sew<2>(cubeLCC.beta(d2, 1), cubeLCC.beta(d8, 1));
	cubeLCC.sew<2>(cubeLCC.beta(d2, 1, 1), cubeLCC.beta(d10, 1));	
	cubeLCC.sew<2>(cubeLCC.beta(d3, 1, 1), cubeLCC.beta(d12, 1, 1));
	cubeLCC.sew<2>(cubeLCC.beta(d3, 1), cubeLCC.beta(d6, 1));
	
	cubeLCC.sew<2>(cubeLCC.beta(d4, 1, 1), cubeLCC.beta(d9, 1, 1));
	cubeLCC.sew<2>(cubeLCC.beta(d5, 1, 1), cubeLCC.beta(d11, 1));
	cubeLCC.sew<2>(cubeLCC.beta(d5, 1), cubeLCC.beta(d7, 1));
	cubeLCC.sew<2>(cubeLCC.beta(d6, 1, 1), cubeLCC.beta(d9, 1));
	cubeLCC.sew<2>(cubeLCC.beta(d7, 1, 1), cubeLCC.beta(d10, 1, 1));
	cubeLCC.sew<2>(cubeLCC.beta(d8, 1, 1), cubeLCC.beta(d11, 1, 1));
	
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

