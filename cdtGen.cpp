#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Delaunay_triangulation_3.h>

using namespace std;


typedef CGAL::Point_3 Point;
typedef CGAL::Delaunay_triangulation_3 Delaunay;

// PLC:
map <unsigned, Point> PLCvertices; // mapping between vertex coordinates and corresponding unique id
list <unsigned, unsigned, unsigned> PLCfaces; // contains ids of vertices/points making the triangle

// Delaunay tetrahedralization:
Delaunay DT;

// CDT: output mesh:
map <unsigned, Point> CDTvertices;
list <unsigned, unsigned, unsigned, unsigned> CDTtets; // contains ids of vertices making tetrahedron


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
	
	list <Point> tempPointList;

	for (map<unsigned, Point>::iterator pit = PLCvertices.begin(); pit != PLCvertices.end(); pit++)
		tempPointList.push(PLCvertices.find(i)->second);

    	DT.insert(tempPointList.begin(). tempPointList.end());
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
