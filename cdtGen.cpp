#include <CGAL/Delaunay_triangulation_3.h>

using namespace std;

// Datastructure implementations:
// PLC: Implement using CGAL Linear cell complexes
// Delaunay tetrahedralization: CGAL::Delaunay_3
// Output mesh: Linear cell complex


void readPLCInput()
{
	// read PLY file(assumed to contain the PLC)
    	// initialize inputVertices
    	// initialize inputFaces


	
}


void computeInitialDelaunayTetrahedralization()
{

    	
}

void removeLocalDegeneracies()
{
	
}

void recoverConstraintFaces()
{
	
}


int main()
{
	readPLCInput();
	computeInitialDelaunayTetrahedralization();
	removeLocalDegeneracies();
	recoverConstraintFaces();
        return 0;
}
