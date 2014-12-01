/*
#include <string.h>
#include <iostream>
#include <unordered_set>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include "rply/rply.h"


using namespace std;
using namespace CGAL;


typedef Exact_predicates_inexact_constructions_kernel K;
typedef Point_3<K> Point;
typedef Delaunay_triangulation_3<K> Delaunay;
typedef Delaunay::Vertex_handle Vertex_handle;


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

/////////////////////////////////////////////// Local Degeneracy Removal ////////////////////////////////////////////////////////////

/*
class DegenerateVertexSetCandidate
{
	Vertex_handle degenSetVertex[5]; // 5 vertices constitute a degeneracy	
	
	bool operator==(const DegenerateVertexSetCandidate &anotherSet)
	{
		bool matched[5];			

		for (unsigned int i = 0; i < 5; i++)
		{	
			matched[i] = false;

			for (unsigned int n = 0; n < 5; n++)
				{
					if (another.degenSetVertex[i] == degenSetVertex[n])
					{
						matched[i] = true;
						break;
					}			
				}
		}

		unsigned int m;
		for (m = 0; m < 5; m++)
		{
			if (matched[m] == true)
				continue;
			else
				break;			
		}

		if (m == 5)
			return true; // duplicate(both key and hash value match)
		else 
			return false; // unique element(only hash values matches)

	}

};

// implements hash function for unorderd_set used to model degeneracy queue and missing face queue
template <>
struct hash<DegenerateVertexSetCandidate>
{
	size_t operator()(const DegenerateVertexSetCandidate &key) const
	{
		size_t hashValue;
		for (unsigned i = 0; i < 5; i++)
			hashValue ^= hash<Vertex_handle>(key.degenVertSet[i]); 

		return (hashValue);
	}
};


*/
	
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

// finds all local degeneracies in DT and adds them to a global queue
void addLocalDegeneraciesToQueue(unordered_set<DegenerateVertexSetCandidate> localDegeneracySet)
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
					{
						if (localDegeneracySet.find(degenerateSet) != localDegeneracySet.end())
							localDegeneracySet.insert(degenerateSet);	
					}
				}
	}
}


/*
// perturbation should be such that it does not make PLC inconsistent
bool isVertexPerturbable(Vertex)
{
	bool pertubable = false;

	
	return perturbable;
}

bool isVertexSegmentSafePertubable()
{
	bool segmentSafePerturbable = false;

	return segmentSafePerturbable;
}

bool isDegeneracyRemovable()
{
	bool removable = false;
	
	if (vertexIsPerturbable(vertex))
		if (vertexIsSegmentSafePerturbable(vertex))
			removable = true;

	return removable;
}
*/

// removes local degeneracies from Delaunay tetrahedralization
void removeLocalDegeneracies()
{
	cout << "\nStarting local degeneracy removal...";

	// compute all local degeneracies in DT and add them to Q
	unordered_set<DegenerateVertexSetCandidate, hashFunction> localDegeneracySet;
	addLocalDegeneraciesToQueue(localDegeneracySet);

	// repeat while Q != NULL	
	while (localDegenercySet.size() != 0)
		for (unordered_set<localDegeneracySet, hashFunction>::iterator qIter = localDegeneracySet.begin(); qIter != localDegeneracySet.end(); qIter++)
			{
				if (isDegeneracyRemovable(qIter))
					perturbRemove(qIter, localDegeneracySet); 
				else
				{
					Vertex vb;	
					computeBreakPoint(vb, qIter, localDegeneracySet);
					if (isEncroachingPLC(vb))
						boundaryProtection(vb);
					else
						inputPLCVertices.push(vb);		
				}						
			}

	cout << "Local degeneracy removal completed";
}

///////////////////////////////////////////////////// Local Degeneracy Removal //////////////////////////////////////////////////////

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
