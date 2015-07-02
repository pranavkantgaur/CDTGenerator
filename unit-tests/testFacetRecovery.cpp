#include "gtest/gtest.h"
#include "cdtGen.h"


/*! \class TestRecoverConstraintFacets
 *  \brief This class provides object for testing implementation of facet recovery. 
 */
class TestRecoverConstraintFacets:public CDTGenerator
{
	public: 
		bool runTest();

	protected:			
		size_t numberOfAlreadyPresentConstraintFacets();
		size_t numberOfRecoveredConstraintFacets();
		bool areFacetsSameGeometry(LCC::Dart_handle, LCC, LCC::Dart_handle, LCC);
};


/*! \fn size_t TestRecoverConstraintFacets::numberOfAlreadyPresentConstraintFacets()
 *  \brief Computes number of constraint facets already present in DT before applying recovery procedure.
 *  \return Number of constraint facets already present in DT.
 */
size_t TestRecoverConstraintFacets::numberOfAlreadyPresentConstraintFacets()
{
	size_t nConstraintFacetsAlreadyPresent = 0;
	Delaunay::Vertex_handle vh1, vh2, vh3;
	Delaunay::Cell_handle c;
	CGALPoint p[3];
	int i, j, k;

	for (LCC::One_dart_per_cell_range<2>::iterator fIter = plc.one_dart_per_cell<2>().begin(), fIterEnd = plc.one_dart_per_cell<2>().end(); fIter != fIterEnd; fIter++)
	{
		for (LCC::One_dart_per_incident_cell_range<0, 2>::iterator pIter = plc.one_dart_per_incident_cell<0, 2>(fIter).begin(), pIterEnd = plc.one_dart_per_incident_cell<0, 2>(fIter).end(); pIter != pIterEnd; pIter++)
			p[i] = plc.point(pIter);

		if (DT.is_vertex(p[0], vh1))
			if (DT.is_vertex(p[1], vh2))
				if (DT.is_vertex(p[2], vh3))
					if (DT.is_facet(vh1, vh2, vh3, c, i, j, k))
						nConstraintFacetsAlreadyPresent++;

	}
	return nConstraintFacetsAlreadyPresent;
}


/*! \fn bool TestRecoverConstraintFacets::areFacetsSameGeometry(LCC::Dart_handle d1, LCC lcc1, LCC::Dart_handle d2, LCC lcc2)
 *  \brief Tests whether input facets from given LCCs are _geometrically_ same.
 *  \param [in] d1 Dart handle to first facet.
 *  \param [in] lcc1 LCC containing first facet.
 *  \param [in] d2 Dart handle to second facet.
 *  \param [in] lcc2 LCC containing second facet.
 */
bool TestRecoverConstraintFacets::areFacetsSameGeometry(LCC::Dart_handle d1, LCC lcc1, LCC::Dart_handle d2, LCC lcc2)
{
	LCC lcc;
	CGALPoint facet1Points[3];
	CGALPoint facet2Points[3];
	size_t i = 0, j = 0;

	for (LCC::One_dart_per_incident_cell_range<0, 2>::iterator pIter1 = lcc1.one_dart_per_incident_cell<0, 2>(d1).begin(), pIter1End = lcc1.one_dart_per_incident_cell<0, 2>(d1).end(); pIter1 != pIter1End; pIter1++)
		facet1Points[i++] = lcc1.point(pIter1);
	
	for (LCC::One_dart_per_incident_cell_range<0, 2>::iterator pIter2 = lcc2.one_dart_per_incident_cell<0, 2>(d2).begin(), pIter2End = lcc2.one_dart_per_incident_cell<0, 2>(d2).end(); pIter2 != pIter2End; pIter2++)
		facet2Points[j++] = lcc2.point(pIter2);

	lcc.make_triangle(facet1Points[0], facet1Points[1], facet1Points[2]);
	lcc.make_triangle(facet2Points[0], facet2Points[1], facet2Points[2]);
	
	lcc.sew3_same_facets();
	
	size_t nFacetsInLCC = 0;
	for (LCC::One_dart_per_cell_range<2>::iterator fIter = lcc.one_dart_per_cell<2>().begin(), fIterEnd = lcc.one_dart_per_cell<2>().end(); fIter != fIterEnd; fIter++)
		nFacetsInLCC++;

	if (nFacetsInLCC == 1)
		return true;
	else
		return false;
}


/*! \fn size_t TestRecoveryConstraintFacets::numberOfRecoveryConstrainedFacets()
 *  \brief Computes number of recovered constraint facets.
 *  \return Number of recovered constraint facets.
 */
size_t TestRecoverConstraintFacets::numberOfRecoveredConstraintFacets()
{
	size_t nConstraintFacetsRecovered = 0;
	
	for (LCC::One_dart_per_cell_range<2>::iterator plcfIter = plc.one_dart_per_cell<2>().begin(), plcfIterEnd = plc.one_dart_per_cell<2>().end(); plcfIter != plcfIterEnd; plcfIter++)
	{
		// search for this facet in original mesh
		for (LCC::One_dart_per_cell_range<2>::iterator cdtMeshfIter = cdtMesh.one_dart_per_cell<2>().begin(), cdtMeshfIterEnd = cdtMesh.one_dart_per_cell<2>().end(); cdtMeshfIter != cdtMeshfIterEnd; cdtMeshfIter++)
			if (areFacetsSameGeometry(plcfIter, plc, cdtMeshfIter, cdtMesh)) 
				nConstraintFacetsRecovered++;
	}	

	return nConstraintFacetsRecovered;
}


/*! \fn bool TestRecoverConstraintFacets::runTest()
 *  \brief Runs the test to check if all constraint facets are recovered.
 *  \return True if all constraint facets are recovered in output mesh.
 */
bool TestRecoverConstraintFacets::runTest()
{
	readPLCInput();
	computeDelaunayTetrahedralization(-1);
	recoverConstraintSegments();
	removeLocalDegeneracies();
	
	size_t nConstraintFacetsAlreadyPresent = 0, nConstraintFacetsRecovered = 0;

	// already present 
	nConstraintFacetsAlreadyPresent = numberOfAlreadyPresentConstraintFacets();	

	recoverConstraintFacets();
	
	// recovered
	nConstraintFacetsRecovered = numberOfRecoveredConstraintFacets();

	// count number of constraint facets recovered in output mesh(LCC) after recovery execution of recovery algorithm.
	cout << "Number of constraint facets already present: " << nConstraintFacetsAlreadyPresent << endl;
	cout << "Number of constraint facets recovered: " << nConstraintFacetsRecovered << endl;
	
	return (nConstraintFacetsAlreadyPresent == nConstraintFacetsRecovered) ? true : false;
}


TEST (ConstraintFacetsRecoveryTest, allFacetsRecovered)
{
	TestRecoverConstraintFacets t;
	ASSERT_EQ(true, t.runTest());
}

