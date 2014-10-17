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

	cout << "Please enter input filename";
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
