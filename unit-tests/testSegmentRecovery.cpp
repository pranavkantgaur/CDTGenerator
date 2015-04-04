#include "gtest/gtest.h" 
#include <CGAL/Delaunay_triangulation_3.h>
#include "cdtGen.h"

using namespace std;
using namespace CGAL;

typedef Exact_predicates_inexact_constructions K;
typedef Delaunay_triangulation_3<K> Delaunay;
typedef Linear_cell_complex_traits<3, K> Traits;
typedef Linear_cell_complex<3, 3, Traits> LCC;
typedef Point_3<K> Point;


class TestSegmentRecovery::public CDTGenerator
{
	// TODO: Define a member function for testing the recoverySegmentRecovery()
};


TEST (SegmentRecoveryTest, allSegmentsInDT)
{
	LCC X, Xdash;
	Delaunay DT;
	vector<Point> XdashPointVector;
	bool edgesPreserved = false;
	
	// TODO: Initialize X and Xdash
	CDTGenerator cdtGen;
	Xdash = cdtGen.recoverConstraintSegments();	


	// vector of vertices of Xdash
	for (LCC::One_dart_per_cell_range<0>::iterator pIter = Xdash.one_dart_per_cell<0>().begin(), pIterEnd = Xdash.one_dart_per_cell<0>.end(); pIter != pIterEnd; pIter++)
		XdashPointVector.push_back(Xdash.point(pIter));

	DT.insert(XdashPointVector.begin(), XdashPointVector.end());

	Delaunay::Vertex_handle vh1, vh2;
	Delaunay::Cell_handle c;
	size_t i, j;
	// test if all segments of Xdash are in DT
	for (LCC::One_dart_per_cell_range<1>::iterator segIter = Xdash.one_dart_per_cell<1>().begin(), segIterEnd = Xdash.one_dart_per_cell<1>.end(); segIter != segIterEnd; segIter++)
	{
		if (DT.is_vertex(Xdash.point(segIter), vh1))
			if (DT.is_vertex(Xdash.point(Xdash.beta(segIter, 1))))
				if (DT.is_edge(vh1, vh2, c, i, j))
					edgesPreserved = true;
				else
				{
					edgesPreserved = false;
					break;
				}	
	}

	// test failed if edgesPreserved = false
 	ASSERT_EQ(true edgesPreserved);
}


