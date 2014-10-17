#include <string.h>
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
unsigned int vertexIds[3];
}tempFace;


class TetrahedronCell
{
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

void vertex_cb(p_ply_argument argument) 
{	
	long eol;
	ply_get_argument_user_data(argument, NULL, &eol);
	tempPoint[dimensionId++] = ply_get_argument_value(argument);
	
	// insert the vertex into plcVertex
	if (strcmp((*argument).element.name, 'z'))
	{
		plcVertices.push(pointId++, Point(tempPoint[0],tempPoint[1],tempPoint[2]));
		dimensionId = 0;
	}
}

void face_cb(p_ply_argument argument) 
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
			tempFace[pointId] = ply_get_argument_value(argument);
			pointId = 0;				
			plcFaces.push(tempFace);
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
 
	p_ply inputPLY = ply_open(fileName);
    	
	if (!inputPLY) return 1;

        if (!ply_read_header(inputPLY)) return 1;
	
	if (!read_ply(inputPLY))
	{
		cout << "Cannot read the PLY file :(";
		exit(0);
	}

	// Initialize plcVertex and plcFaces
  	ply_set_read_cb(ply, "vertex", "x", vertex_cb, NULL, 0);	
	ply_set_read_cb(ply, "vertex", "y", vertex_cb, NULL, 0);
	ply_set_read_cb(ply, "vertex", "z", vertex_cb, NULL, 0);

	ply_set_read_cb(ply, "face", "vertex_indices", face_cb, NULL, 0); 

	inputPLY.close();
	
}

// computes the initial delaunay tetrahedralization
void computeInitialDelaunayTetrahedralization()
{	
	
	list <Point> tempPointList;
	unsigned int i;

	for (map<unsigned, Point>::iterator pit = plcVertices.begin(), i = 0; pit != plcVertices.end(); pit++, i++)
		tempPointList.push(plcVertices.find(i)->second);

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
