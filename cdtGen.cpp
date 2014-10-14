#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Delaunay_triangulation_3.h>




// reads input PLC
void readPLCInput()
{
	// read PLY file(assumed to contain the PLC)
    	// initialize inputVertices
    	// initialize inputFaces


	
}

// computes the initial delaunay tetrahedralization
void computeInitialDelaunayTetrahedralization()
{

    	
}

// removes local degeneracies from Delaunay tetrahedralization
void removeLocalDegeneracies()
{
	// compute all local degeneracies in DT and add them to Q
	// repeat while Q != NULL
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
	// end repeat
}	

// recovers the constraint faces
void recoverConstraintFaces()
{
	
}

// main procedure
int main()
{
	readPLCInput();
	computeInitialDelaunayTetrahedralization();
	removeLocalDegeneracies();
	recoverConstraintFaces();
        return 0;
}
