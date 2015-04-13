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

	for (Delaunay::Finite_cells_iterator cIter = DT.finite_cells_begin(), cIterEnd = DT.finite_cells_end(); cIter != cIterEnd; cIter++)
	{
		// TODO: Count number of local degeneracies here.	
	}
}


/*! \fn bool TestLocalDegeneracyRemoval :: runTest()
 *  \brief Tests local degeneracy removal code.
 *  \return True if there are no local degeneracies left after _removal_ operation else returns False.
 */
bool TestLocalDegeneracyRemoval :: runTest()
{
	TestSegmentRecovery setTest;
	size_t nLocalDegeneraciesBeforeRemoval = 0, nLocalDegeneraciesAfterRemoval = 0;

	if (segTest.runTest())
	{
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
	else
		return false;
}

TEST (LocalDegeneracyRemovalTest, allDegeneraciesRemoved)
{
	TestLocalDegeneracyRemoval t;
	ASSERT_EQ(true, t.runTest());
}







