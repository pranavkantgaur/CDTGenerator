#include "gtest/gtest.h"
#include "cdtGen.h"

/*! \class TestInitialDT
 *  \brief Provides object for testing implementation of initial Delaunay triangulation step in CDT algorithm.
 */
class TestInitialDT : public CDTGenerator
{
	public:
		bool runTest();
	protected:
		size_t runInitialDT();
};
