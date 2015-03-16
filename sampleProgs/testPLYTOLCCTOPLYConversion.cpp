/*
This code reads a PLY file, generates corresponding LCC, writes LCC back to PLY format. It uses RPly library for interfacing with PLY format.
*/

#include "../rply/rply.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>

using namespace std;
using namespace CGAL;

struct MyItem
{
template<class Refs>
struct Dart_wrapper
{
typedef CGAL::Dart<3, Refs> Dart;
typedef CGAL::Cell_attribute_with_point<Refs, unsigned int, Tag_true> Vertex_attribute;
typedef CGAL::cpp11::tuple<Vertex_attribute> Attributes;
};
};

typedef Exact_predicates_inexact_constructions_kernel K;
typedef Point_3<K> CGALPoint;
typedef Linear_cell_complex_traits<3, K> Traits;
typedef Linear_cell_complex<3, 3, Traits, MyItem> LCC;
typedef LCC::Dart_handle DartHandle;


class Point
{
	public:
		float  x, y, z;
		Point(float aX, float aY, float aZ)
		{
			x = aX;
			y = aY;
			z = aZ;
		}
};

class Triangle
{
	public:
		size_t pointIds[3];
		Triangle(size_t point1Id, size_t point2Id, size_t point3Id)
		{
			pointIds[0] = point1Id;
			pointIds[1] = point2Id;
			pointIds[2] = point3Id;
		}
};

vector<Point> vertexVector;
vector<Triangle> faceVector;
LCC lcc;


static float tempPoint[3];
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
		vertexVector.push_back(Point(tempPoint[0],tempPoint[1],tempPoint[2]));
		dimensionId = 0;
	}
	
	return 1;
}

static int face_cb(p_ply_argument argument) 
{
	long length, value_index;
        static unsigned int tempFacePointIds[3];
	
	ply_get_argument_property(argument, NULL, &length, &value_index);

        switch (value_index) 
	{
        	case 0:
	        case 1: 
        		tempFacePointIds[pointId++] = ply_get_argument_value(argument);
		        break;
        	case 2:	
			tempFacePointIds[pointId] = ply_get_argument_value(argument);
			pointId = 0;				
			faceVector.push_back(Triangle(tempFacePointIds[0], tempFacePointIds[1], tempFacePointIds[2]));
			break;
       		default: 
                	break;
        } 
    
	return 1;
}

// reads input PLC
void readPLYDataset()
{
	// read PLY file(assumed to contain the PLC)
	string fileName;

	cout << "\nPlease enter input filename:\t";
	cin >> fileName;
	
	p_ply inputPLY = ply_open(fileName.c_str(), NULL, 0, NULL);
    	
	if (!inputPLY) exit(0);

        if (!ply_read_header(inputPLY)) exit(0);

	// Initialize plcVertex and plcFaces
  	ply_set_read_cb(inputPLY, "vertex", "x", vertex_cb, NULL, 0);	
	ply_set_read_cb(inputPLY, "vertex", "y", vertex_cb, NULL, 0);
	ply_set_read_cb(inputPLY, "vertex", "z", vertex_cb, NULL, 1);

	ply_set_read_cb(inputPLY, "face", "vertex_indices", face_cb, NULL, 0); 
	if (!ply_read(inputPLY))
	{
		cout << "Cannot read the PLY file :(";
		exit(0);
	}

	ply_close(inputPLY);

	cout << "Number of vertices:" << vertexVector.size() << "\n";
	cout << "Number of faces:" << faceVector.size() << "\n";
}


bool areGeometricallySame(DartHandle d1, DartHandle d2)
{
	
	if (lcc.point(d1) == lcc.point(lcc.beta(d2, 1)))
			if (lcc.point(lcc.beta(d1, 1)) == lcc.point(d2))
				return true;
	return false;
}



void generateLCC()
{
	unsigned int vertexIds[3];
	CGALPoint trianglePoints[3];
	vector <pair<DartHandle, DartHandle> > twoCellsToBeSewed;
	
	// reserve a mark
	int sewedMark = lcc.get_new_mark();

	if (sewedMark == -1)
	{
		cout << "\nNo free mark available!!";
		exit(0);
	}

	for (unsigned int n = 0, m = faceVector.size(); n < m; n++)
	{
		for (unsigned int k = 0; k < 3; k++)
			vertexIds[k] = faceVector[n].pointIds[k];
	
		for (unsigned int i = 0; i < 3; i++)
			trianglePoints[i] = CGALPoint(vertexVector[vertexIds[i]].x, vertexVector[vertexIds[i]].y, vertexVector[vertexIds[i]].z);
			
		lcc.make_triangle(trianglePoints[0], trianglePoints[1], trianglePoints[2]);
	}

	// sew facets sharing edge
	for (LCC::Dart_range::iterator segIter1 = lcc.darts().begin(), segIterEnd1 = lcc.darts().end(); segIter1 != segIterEnd1; segIter1++)
	{
		if (!lcc.is_marked(segIter1, sewedMark)) // not sewed till now
		{
			for (LCC::Dart_range::iterator segIter2 = lcc.darts().begin(), segIterEnd2 = lcc.darts().end(); segIter2 != segIterEnd2; segIter2++)
			{
				if (!lcc.is_marked(segIter2, sewedMark) && lcc.is_sewable<2>(segIter1, segIter2))
				{
					if (areGeometricallySame(segIter1, segIter2)) // checks the geometry of segments
					{
						lcc.mark(segIter1, sewedMark);
						lcc.mark(segIter2, sewedMark);
						twoCellsToBeSewed.push_back(pair<DartHandle, DartHandle>(segIter1, segIter2));
						break;
					}
				}
				else
					continue;
			}
		}	
		else
			continue;
	}
	
	// sew the faces sharing an edge
	unsigned int k = 0;
	for (vector<pair<DartHandle, DartHandle> >::iterator dIter = twoCellsToBeSewed.begin(), dIterEnd = twoCellsToBeSewed.end(); dIter != dIterEnd; dIter++)
		if (lcc.is_sewable<2>(dIter->first, dIter->second))
		{
			lcc.sew<2>(dIter->first, dIter->second);
			k++;	
		}

	cout << "\nNumber of sewable facets: " << k;
	cout << "\nNumber of pairs to be sewed: " << twoCellsToBeSewed.size() << "\n";

	return;
}

void testLCC()
{
	unsigned int nVertices = 0;
	for (LCC::One_dart_per_cell_range<0>::iterator vIter = lcc.one_dart_per_cell<0>().begin(); vIter != lcc.one_dart_per_cell<0>().end(); vIter++)
		nVertices++;

	cout << "\nNumber of vertices: " << nVertices;	

	
	unsigned int nSegments = 0;
	for (LCC::One_dart_per_cell_range<1>::iterator sIter = lcc.one_dart_per_cell<1>().begin(); sIter != lcc.one_dart_per_cell<1>().end(); sIter++)
		nSegments++;

	cout << "\nNumber of segments: " << nSegments;	

	unsigned int nFaces = 0;
	for (LCC::One_dart_per_cell_range<2>::iterator fIter = lcc.one_dart_per_cell<2>().begin(); fIter != lcc.one_dart_per_cell<2>().end(); fIter++)
		nFaces++;

	cout << "\nNumber of faces: " << nFaces << "\n";	
}

void writePLYOutput()
{

	p_ply lccOutputPLY;

	if ((lccOutputPLY = ply_create("lccOutputPLY.ply", PLY_ASCII, NULL, 0, NULL)) == NULL)
	{
		cout << "\nCannot open file for writing!!";
		exit(0);
	}

	// count number of vertices and faces in LCC // TODO: No function in Linear cell complex
	size_t nVertices = 0, nFaces = 0;
	for (LCC::One_dart_per_cell_range<0>::iterator pointCountIter = lcc.one_dart_per_cell<0>().begin(), pointCountIterEnd = lcc.one_dart_per_cell<0>().end(); pointCountIter != pointCountIterEnd; pointCountIter++)
		nVertices++;

	for (LCC::One_dart_per_cell_range<2>::iterator faceCountIter = lcc.one_dart_per_cell<2>().begin(), faceCountIterEnd = lcc.one_dart_per_cell<2>().end(); faceCountIter != faceCountIterEnd; faceCountIter++)
		nFaces++;

	ply_add_element(lccOutputPLY, "vertex", nVertices);
	ply_add_scalar_property(lccOutputPLY, "x", PLY_FLOAT);
	ply_add_scalar_property(lccOutputPLY, "y", PLY_FLOAT);
	ply_add_scalar_property(lccOutputPLY, "z", PLY_FLOAT);

	ply_add_element(lccOutputPLY, "face", nFaces);
	ply_add_list_property(lccOutputPLY, "vertex_indices", PLY_UCHAR, PLY_INT32);

	if (!ply_write_header(lccOutputPLY))
	{
		cout << "Header cannot be writen!!";
		exit(0);
	}

	// write vertices
	size_t pointId = 0;
	for (LCC::One_dart_per_cell_range<0>::iterator pointIter = lcc.one_dart_per_cell<0>().begin(), pointIterEnd = lcc.one_dart_per_cell<0>().end(); pointIter != pointIterEnd; pointIter++)
	{
		CGALPoint pt = lcc.point(pointIter); 
		
		ply_write(lccOutputPLY, pt.x());
		ply_write(lccOutputPLY, pt.y());
		ply_write(lccOutputPLY, pt.z());

		lcc.info<0>(pointIter) = pointId++;
	}
	
        // write polygons	
	for (LCC::One_dart_per_cell_range<2>::iterator faceIter = lcc.one_dart_per_cell<2>().begin(), faceIterEnd = lcc.one_dart_per_cell<2>().end(); faceIter != faceIterEnd; faceIter++)
	{
		ply_write(lccOutputPLY, 3);
		for (LCC::One_dart_per_incident_cell_range<0, 2>::iterator pointInFaceIter = lcc.one_dart_per_incident_cell<0, 2>(faceIter).begin(), pointInFaceIterEnd = lcc.one_dart_per_incident_cell<0, 2>(faceIter).end(); pointInFaceIter != pointInFaceIterEnd; pointInFaceIter++)
			ply_write(lccOutputPLY, lcc.info<0>(pointInFaceIter)); // TODO: How to ensure that it is writen in counter clockwise order only?		 
	}

	ply_close(lccOutputPLY);

}



int main()
{
	readPLYDataset();
	generateLCC();
	testLCC();
	writePLYOutput();
	return 0;
}
