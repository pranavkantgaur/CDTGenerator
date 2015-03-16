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

LCC cubeLCC;



bool areGeometricallySame(DartHandle d1, DartHandle d2)
{
	
	if (cubeLCC.point(d1) == cubeLCC.point(cubeLCC.beta(d2, 1)))
			if (cubeLCC.point(cubeLCC.beta(d1, 1)) == cubeLCC.point(d2))
				return true;
	return false;
}



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

	// add triangles representing cube to cubeLCC
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

	vector <pair<DartHandle, DartHandle> > twoCellsToBeSewed;
	int sewedMark = cubeLCC.get_new_mark();	

	// sew facets sharing edge
	for (LCC::Dart_range::iterator segIter1 = cubeLCC.darts().begin(), segIterEnd1 = cubeLCC.darts().end(); segIter1 != segIterEnd1; segIter1++)
	{
		if (!cubeLCC.is_marked(segIter1, sewedMark)) // not sewed till now
			{
					
				for (LCC::Dart_range::iterator segIter2 = cubeLCC.darts().begin(), segIterEnd2 = cubeLCC.darts().end(); segIter2 != segIterEnd2; segIter2++)
				{
					if (!cubeLCC.is_marked(segIter2, sewedMark) && cubeLCC.is_sewable<2>(segIter1, segIter2))
					{
						if (areGeometricallySame(segIter1, segIter2)) // checks the geometry of segments
						{
							cubeLCC.mark(segIter1, sewedMark);
							cubeLCC.mark(segIter2, sewedMark);
							twoCellsToBeSewed.push_back(pair<DartHandle, DartHandle>(segIter1, segIter2));
							break;
						}
					}
					else
						continue;
				}
			}	
		else
			continue;
	}

	// sew the faces sharing an edge
	unsigned int k = 0;
	for (vector<pair<DartHandle, DartHandle> >::iterator dIter = twoCellsToBeSewed.begin(), dIterEnd = twoCellsToBeSewed.end(); dIter != dIterEnd; dIter++)
		if (cubeLCC.is_sewable<2>(dIter->first, dIter->second))
		{
			cubeLCC.sew<2>(dIter->first, dIter->second);
			k++;	
		}

	// Verify correctness
	cout << "\nNumber of sewable facets: " << k;
	cout << "\nNumber of pairs to be sewed: " << twoCellsToBeSewed.size() << "\n";

	unsigned int nVertices = 0;
	for (LCC::One_dart_per_cell_range<0>::iterator vIter = cubeLCC.one_dart_per_cell<0>().begin(); vIter != cubeLCC.one_dart_per_cell<0>().end(); vIter++)
		nVertices++;

	cout << "\nNumber of vertices: " << nVertices;	

	
	unsigned int nSegments = 0;
	for (LCC::One_dart_per_cell_range<1>::iterator sIter = cubeLCC.one_dart_per_cell<1>().begin(); sIter != cubeLCC.one_dart_per_cell<1>().end(); sIter++)
		nSegments++;

	cout << "\nNumber of segments: " << nSegments;	

	unsigned int nFaces = 0;
	for (LCC::One_dart_per_cell_range<2>::iterator fIter = cubeLCC.one_dart_per_cell<2>().begin(); fIter != cubeLCC.one_dart_per_cell<2>().end(); fIter++)
		nFaces++;

	cout << "\nNumber of faces: " << nFaces << "\n";	


	cubeLCC.free_mark(sewedMark); 


	return 0;
	
}

