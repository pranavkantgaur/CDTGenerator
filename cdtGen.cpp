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
	// compute all local degeneracies in DT and add them to Q
	// for each local degeneracy l in Q:
		// if l is removable:
			// remove l by small perturbation
		// else,
			// compute vb to break l
				// if vb encroaches upon any segment or subface,
					// push l into queue
					// call boundary protection procedure
				// else
					// insert vb to break l 	
				// end if
		// end if
	// end for
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
