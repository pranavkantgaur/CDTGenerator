#include "gtest/gtest.h"
#include "cdtGen.h"

/*! \class TestLocalDegeneracyRemoval
*   \brief Provides object for testing implementation of local degeneracy removal process. 
*/
class TestLocalDegeneracyRemoval : public CDTGenerator
{
	public:
		bool runTest();
	protected:
		size_t locateLocalDegeneracies();	
};

/*! \fn TestLocalDegeneracyRemoval :: locateLocalDegeneracy()
 *  \brief Computes number of local degeneracies in PLC.
 *
 *  Computes local degeneracies in PLC by counting number of _cospherical_ 5-point sets in corresponding Delaunay triangulation.
 *  \return Returns number of local degenracies in PLC.
 */
size_t TestLocalDegeneracyRemoval :: locateLocalDegeneracies()
{
	size_t nLocalDegeneracies = 0;
	CGALPoint p[5];

	for (Delaunay::Finite_cells_iterator cIter = DT.finite_cells_begin(), cIterEnd = DT.finite_cells_end(); cIter != cIterEnd; cIter++)
	{
		for (size_t i = 0; i < 4; i++)
			p[i] = cIter->vertex(i)->point();			
		for (size_t neighborID = 0; neighborID < 4; neighborID++)
		{
			Delaunay::Cell_handle neighborCell = cIter->neighbor(neighborID);
			p[4] = (DT.mirror_vertex(cIter, neighborID))->point();
			if (side_of_bounded_sphere(p[0], p[1], p[2], p[3], p[4]) == ON_BOUNDARY)
				nLocalDegeneracies++;
		}
	}
	if (nLocalDegeneracies % 2 == 0)
		cout << "OK";
	else
		cout << "ERROR";

	return nLocalDegeneracies;
} 


/*! \fn bool TestLocalDegeneracyRemoval :: runTest()
 *  \brief Tests local degeneracy removal code.
 *  \return True if there are no local degeneracies left after _removal_ operation else returns False.
 */
bool TestLocalDegeneracyRemoval :: runTest()
{
	readPLCInput();
	computeDelaunayTetrahedralization();
	recoverConstraintSegments();

	size_t nLocalDegeneraciesBeforeRemoval = 0, nLocalDegeneraciesAfterRemoval = 0;

	// Now we can check for correctness of local degeneracy removal
	nLocalDegeneraciesBeforeRemoval = locateLocalDegeneracies();
	removeLocalDegeneracies(); // plc is modified
	nLocalDegeneraciesAfterRemoval = locateLocalDegeneracies();	
	cout << "Number of local Degeneracies before removal: " << nLocalDegeneraciesBeforeRemoval << endl;
	cout << "Number of local Degeneracies after removal: " << nLocalDegeneraciesAfterRemoval << endl;

	if (nLocalDegeneraciesAfterRemoval != 0)
		return false;
	else 
		return true;	
}

TEST (LocalDegeneracyRemovalTest, allDegeneraciesRemoved)
{
	TestLocalDegeneracyRemoval t;
	ASSERT_EQ(true, t.runTest());
}







