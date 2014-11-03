/*
#include <string.h>
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include "rply/rply.h"


using namespace std;
using namespace CGAL;


typedef Exact_predicates_inexact_constructions_kernel K;
typedef Point_3<K> Point;
typedef Delaunay_triangulation_3<K> Delaunay;

// PLC:
map <unsigned, Point> plcVertices; // mapping between vertex coordinates and corresponding unique id

class TriangleFace
{
	public:
		unsigned int vertexIds[3];
}tempFace;


class TetrahedronCell
{
	public:
		unsigned int vertexIds[4];
}tempTet;

list <TriangleFace> plcFaces; // contains ids of vertices/points making the triangle

// Delaunay tetrahedralization:
Delaunay DT;

// CDT: 
map <unsigned, Point> cdtVertices;
list <TriangleFace> cdtTets; // contains ids of vertices making tetrahedron

static unsigned tempPoint[3];
unsigned int dimensionId = 0;
static unsigned pointId = 0;

static int vertex_cb(p_ply_argument argument) 
{	
	long eol;
	ply_get_argument_user_data(argument, NULL, &eol);
	tempPoint[dimensionId++] = ply_get_argument_value(argument);
	
	// insert the vertex into plcVertex
	if (eol)
	{
		plcVertices[pointId++] = Point(tempPoint[0],tempPoint[1],tempPoint[2]);
		dimensionId = 0;
	}
}

static int face_cb(p_ply_argument argument) 
{
	long length, value_index;

	ply_get_argument_property(argument, NULL, &length, &value_index);

        switch (value_index) 
	{
        	case 0:
	        case 1: 
        		tempFace.vertexIds[pointId++] = ply_get_argument_value(argument);
	                break;
        	case 2:	
			tempFace.vertexIds[pointId] = ply_get_argument_value(argument);
			pointId = 0;				
			plcFaces.push_front(tempFace);
	                break;
        	default: 
                	cout << "Invalid number of points in a face specified :(";
			break;
        } 
    
}

// reads input PLC
void readPLCInput()
{
	// read PLY file(assumed to contain the PLC)
	string fileName;

	cout << "\nPlease enter input filename:\t";
	cin >> fileName;
	
	p_ply inputPLY = ply_open(fileName.c_str(), NULL, 0, NULL);
    	
	if (!inputPLY) exit(0);

        if (!ply_read_header(inputPLY)) exit(0);
	
	if (!ply_read(inputPLY))
	{
		cout << "Cannot read the PLY file :(";
		exit(0);
	}

	// Initialize plcVertex and plcFaces
  	ply_set_read_cb(inputPLY, "vertex", "x", vertex_cb, NULL, 0);	
	ply_set_read_cb(inputPLY, "vertex", "y", vertex_cb, NULL, 0);
	ply_set_read_cb(inputPLY, "vertex", "z", vertex_cb, NULL, 1);

	ply_set_read_cb(inputPLY, "face", "vertex_indices", face_cb, NULL, 0); 

	ply_close(inputPLY);
	
}

// computes the initial delaunay tetrahedralization
void computeInitialDelaunayTetrahedralization()
{	
	
	list <Point> tempPointList;
	unsigned int i = 0;

	for (map<unsigned, Point>::iterator pit = plcVertices.begin(); pit != plcVertices.end(); pit++, i++)
		tempPointList.push_front(plcVertices.find(i)->second);

    	DT.insert(tempPointList.begin(), tempPointList.end());

	cout << "\nInitial Delaunay tetrahedralization computed!!";
}
*/

class DegenerateVertexSetCandidate
{
	Vertex_handle degenSetVertex[5];	
};


// returns true if  given vertices are co-spherical
bool areCospherical(DegenerateVertexSetCandidate degenSet)
{
	
	Point_3 p[5];
	
	for (unsigned int i = 0; i < 5; i++)	
		p[i] = (degenSet.degenSetVertex[i])->point();

	if (CGAL::side_of_bounded_sphere(p[0],p[1],p[2],p[3],p[4]) == CGAL::ON_BOUNDARY)
		return true;
	else
		return false;
}

// removes duplicate vertex sets from global degeneracyQueue
void removeDuplicateDegenracies()
{

}


// finds all local degeneracies in DT and adds them to a global queue
void addLocalDegeneraciesToQueue(queue<DegenerateVertexSetCandidate> degeneracyQueue)
{
	
	DegenerateVertexSetCandidate degenerateSet;
	
	for (delaunay::Finite_cell_iterator cellIter = DT.finite_cells_begin(); cellIter != DT.finite_cells_end(); cellIter++)
	{
		for (unsigned int n = 0; n < 4; n++)
			degenerateSet.degenSetVertex[n] = cellIter->vertex(n);

		for (unsigned int j = 0; j < 4; j++)
			for (unsigned int k = 0; k < 4; k++)
				{
					if ((cellId->neighbor(j))->neighbor(k) == cellId)
						degenerateSet.degenSetVertex[4] = (cellIter->neighbor)->vertex(k);		
						
					if (areCospherical(degenerateSet))
						degeneracyQueue.push_front(degenerateSet);	
				}
	}

	removeDuplicateDegeneracies(); 
}



// removes local degeneracies from Delaunay tetrahedralization
void removeLocalDegeneracies()
{
	
	cout << "\nStarting local degeneracy removal...";

	// compute all local degeneracies in DT and add them to Q
	queue<DegenerateVertexSetCandidate> degeneracyQueue;
	addLocalDegeneraciesToQueue(degeneracyQueue);

	// repeat while Q != NULL	
	while (degenercyQueue.size() != 0)
		for (queue<>::iterator qIter = degeneracyQueue.begin(); qIter != degeneracyQueue.end(); qIter++)
			{
				if (isDegeneracyRemovable(qIter))
					perturbRemove(qIter, degeneracyQueue); 
				else
				{
					Vertex vb;	
					computeBreakPoint(vb, qIter, degeneracyQueue);
					if (isEncroachingPLC(vb))
						boundaryProtection(vb);
					else
						inputPLCVertices.push(vb);		
				}						
			}

	cout << "Local degeneracy removal completed";
}


/*
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
*/
