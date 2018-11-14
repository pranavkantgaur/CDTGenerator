#include "gtest/gtest.h"
#include "cdtGen.h"

/*! \class TestPLCInput
 *  \brief Derived from CDTGenerator class, this class defines unit-tests for evaluating the logic of loading input file and (internally)representing it as Linear Cell Complex.
 */
class TestPLCInput:public CDTGenerator
{
        public:
        bool runTest();          
}
