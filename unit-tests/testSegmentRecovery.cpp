#include "gtest/gtest.h" 
#include "cdtGen.h"


/*! \class TestSegmentRecovery
 *  \brief Derived from CDTGenerator class, this class defines a unit-test for segment recovery.
 */
class TestSegmentRecovery:public CDTGenerator
{
	public:
		bool runTest();	
};

/*! \fn bool TestSegmentRecovery::runTest()
*   \brief Tests segment recovery function for CDTGenerator class.
*   \return ```True``` if all constraint subsegments are recovered otherwise returns ```False```.
*   Tests whether all segments of original PLC are present in final PLC even as union of subsegments. It also tests whether all constriant subsegments of final PLC are recovered in final Delaunay tetrahedralization ```DT```.
*   
*/
bool TestSegmentRecovery::runTest()
{
	readPLCInput();
	computeDelaunayTetrahedralization();

	Delaunay::Vertex_handle vh1, vh2;
	Delaunay::Cell_handle c;
	int i, j;
	size_t nEdgesAlreadyPresent = 0;

	for (LCC::One_dart_per_cell_range<1>::iterator segIter = plc.one_dart_per_cell<1>().begin(), segIterEnd = plc.one_dart_per_cell<1>().end(); segIter != segIterEnd; segIter++)
	{
		if (DT.is_vertex(plc.point(segIter), vh1))
			if (DT.is_vertex(plc.point(plc.beta(segIter, 1)), vh2))
				if (DT.is_edge(vh1, vh2, c, i, j))
					nEdgesAlreadyPresent++;
	}

	recoverConstraintSegments();
	
	vector<CGALPoint> plcVertexVector;
	
	// vector of vertices of plc
	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
		plcVertexVector.push_back(plc.point(pIter));

	DT.insert(plcVertexVector.begin(), plcVertexVector.end());


	// test if all segments of plc are in DT after executing segment recovery procedure
	size_t nEdgesPresentAfterSegRecovery = 0, nEdgesNotPresentEvenAfterSegRecovery = 0;

	for (LCC::One_dart_per_cell_range<1>::iterator segIter = plc.one_dart_per_cell<1>().begin(), segIterEnd = plc.one_dart_per_cell<1>().end(); segIter != segIterEnd; segIter++)
	{
		if (DT.is_vertex(plc.point(segIter), vh1))
			if (DT.is_vertex(plc.point(plc.beta(segIter, 1)), vh2))
				if (DT.is_edge(vh1, vh2, c, i, j))
					nEdgesPresentAfterSegRecovery++;				
				else
					nEdgesNotPresentEvenAfterSegRecovery++;
		
	}

	cout << "Number of edges already present: " << nEdgesAlreadyPresent << "\n";
	cout << "Number of edges present after segment recovery: " << nEdgesPresentAfterSegRecovery << "\n";
	cout << "Number of edges not present even after segment recovery: " << nEdgesNotPresentEvenAfterSegRecovery << "\n"; 


	return (nEdgesNotPresentEvenAfterSegRecovery == 0) ? true : false;
}


TEST (SegmentRecoveryTest, allSegmentsInDT)
{
	TestSegmentRecovery t;
 	ASSERT_EQ(true, t.runTest());
}


