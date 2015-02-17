/* 
 * Minimal code using CGAL's LCC package for constructing linear cell complex of a cube 
 */
#include <vector>
#include <tuple>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Linear_cell_complex.h>

using namespace std;
using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;

int main()
{
	vector<tuple<unsigned int, unsigned int, unsigned int> > cubePoints;

	// create storage for points of a cube
	cubePoints.push_back(make_tuple(0, 0, 0)); // points of a cube
	cubePoints.push_back(make_tuple(0, 0, 1));
	cubePoints.push_back(make_tuple(0, 1, 0));
	cubePoints.push_back(make_tuple(0, 1, 1));
	cubePoints.push_back(make_tuple(1, 0, 0));
	cubePoints.push_back(make_tuple(1, 0, 1));
	cubePoints.push_back(make_tuple(1, 1, 0));
	cubePoints.push_back(make_tuple(1, 1, 1));

	LCC cubeLCC;

	// add triangles representing cube to lcc
	cubeLCC.make_triangle(cubePoints[]);

	// sew triangles with identical geometry
	
	// Verify correctness
	 
	return 0;
}

