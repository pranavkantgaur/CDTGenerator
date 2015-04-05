#include "gtest/gtest.h" 
#include "cdtGen.h"

using namespace std;
using namespace CGAL;

typedef Exact_predicates_inexact_constructions K;
typedef Delaunay_triangulation_3<K> Delaunay;
typedef Linear_cell_complex_traits<3, K> Traits;
typedef Linear_cell_complex<3, 3, Traits> LCC;

/*! \class TestSegmentRecovery
 *  \brief Derived from CDTGenerator class, this class defines a unit-test for segment recovery.
 */
class TestSegmentRecovery::public CDTGenerator
{
	public:
		bool runTest();	
	
	protected:
		bool edgesPreserved;
};


bool TestSegmentRecovery::runTest()
{
	readInput();
	recoverConstraintSegments();
	
	plcVertexVector.clear();
	
	// vector of vertices of plc
	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>.end(); pIter != pIterEnd; pIter++)
		plcVertexVector.push_back(plc.point(pIter));

	DT.insert(plcVertexVector.begin(), plcVertexVector.end());

	Delaunay::Vertex_handle vh1, vh2;
	Delaunay::Cell_handle c;
	size_t i, j;
	// test if all segments of plc are in DT
	for (LCC::One_dart_per_cell_range<1>::iterator segIter = plc.one_dart_per_cell<1>().begin(), segIterEnd = plc.one_dart_per_cell<1>.end(); segIter != segIterEnd; segIter++)
	{
		if (DT.is_vertex(plc.point(segIter), vh1))
			if (DT.is_vertex(plc.point(plc.beta(segIter, 1))))
				if (DT.is_edge(vh1, vh2, c, i, j))
					edgesPreserved = true;
				else
				{
					edgesPreserved = false;
					break;
				}	
	}
	return edgesPreserved;
}


TEST (SegmentRecoveryTest, allSegmentsInDT)
{
	TestSegmentRecovery t;
 	ASSERT_EQ(true t.runTest());
}


