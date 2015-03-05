#include <iostream>
#include <vector>
#include <math.h>
#include <map>
#include <unordered_set>

#include <CGAL/Random.h>
#include <CGAL/Object.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Circle_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Spherical_kernel_intersections.h>
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/enum.h>

#include "rply/rply.h"


#define Pi 22.0/7.0
#define INVALID_VALUE -1.0f // used in context of distances 
using namespace std;
using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;
typedef Linear_cell_complex_traits<3, K> Traits;
typedef Linear_cell_complex<3, 3, Traits> LCC;
typedef LCC::Dart_handle DartHandle;
typedef Triangulation_vertex_base_with_info_3<DartHandle, K> Vb; 
typedef Triangulation_data_structure_3<Vb> Tds;
typedef Delaunay_triangulation_3<K, Tds, Fast_location> Delaunay;
typedef Delaunay::VertexHandle VertexHandle;
typedef Delaunay::Cell_handle CellHandle;
typedef Point_3<K> CGALPoint;
typedef Circle_3<K> CGALCircle;
typedef Sphere_3<K> CGALSphere;
typedef Tetrahedron_3<K> CGALTetrahedron;
typedef Triangle_3<K> CGALTriangle;

typedef Exact_spherical_kernel_3 SK;
typedef Line_arc_3<SK> CGALSphericalSegment;
typedef Sphere_3<SK> CGALSphericalSphere;
typedef Segment_3<SK> CGALSphericalSegment; 
typedef Point_3<SK> CGALSphericalPoint;



/*
 * Input  : PLC(represented using CGAL's LCC structure)
 * Output : cdtMesh(again, LCC)
 */

// Input
LCC plc;

// Intermidiate global structures
Delaunay DT;
vector <Point> plcVertexVector;
class Triangle
{
	size_t pointIds[3];
};

vector <Triangle> plcFaceVector;

// Output
LCC cdtMesh;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


static float tempPoint[3];
static int pointCount = 0;
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
		plcVertexVector.push_back(Point(tempPoint[0],tempPoint[1],tempPoint[2]));
		dimensionId = 0;
	}
	
	return 1;
}

static int face_cb(p_ply_argument argument) 
{
	long length, value_index;
        static Triangle tempFace;
	
	ply_get_argument_property(argument, NULL, &length, &value_index);

        switch (value_index) 
	{
        	case 0:
	        case 1: 
        		tempFace.pointIds[pointId++] = ply_get_argument_value(argument);
		        break;
        	case 2:	
			tempFace.pointIds[pointId] = ply_get_argument_value(argument);
			pointId = 0;				
			plcFaceVector.push_back(tempFace);
			break;
       		default: 
                	break;
        } 
    
	return 1;
}



bool areGeometricallySameSegments(DartHandle d1, DartHandle d2, LCC &lcc)
{
	if (lcc.point(d1) == lcc.point(lcc.beta(d2, 1))
			if (lcc.point(lcc.beta(d1, 1)) == lcc.point(d2))
				return true;

	return false;				
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

	cout << "Number of vertices:" << plcVertexVector.size() << "\n";
	cout << "Number of faces:" << plcFaceVector.size() << "\n";

	// Initialize PLC
	vector<pair<DartHandle, DartHandle> > twoCellsToBeSewed;
	size_t vertexIds[3];
 
	int sewedMark = lcc.get_new_mark();

	if (sewedMark == -1)
	{
		cout << "\nNo free mark available!!";
		exit(0);
	}

	for (unsigned int n = 0, m = faceVector.size(); n < m; n++)
	{
		for (unsigned int k = 0; k < 3; k++)
			vertexIds[k] = plcFaceVector[n].pointIds[k];
	
		for (unsigned int i = 0; i < 3; i++)
			trianglePoints[i] = CGALPoint(plcVertexVector[vertexIds[i]].x(), plcVertexVector[vertexIds[i]].y(), plcVertexVector[vertexIds[i]].z());
			
		plc.make_triangle(trianglePoints[0], trianglePoints[1], trianglePoints[2]);
	}

	// sew facets sharing edge
	for (LCC::Dart_range::iterator segIter1 = plc.darts().begin(), segIterEnd1 = plc.darts().end(); segIter1 != segIterEnd1; segIter1++)
	{
		if (!plc.is_marked(segIter1, sewedMark)) // not sewed till now
		{
			for (LCC::Dart_range::iterator segIter2 = plc.darts().begin(), segIterEnd2 = plc.darts().end(); segIter2 != segIterEnd2; segIter2++)
			{
				if (!plc.is_marked(segIter2, sewedMark) && plc.is_sewable<2>(segIter1, segIter2))
				{
					if (areGeometricallySameSegments(segIter1, segIter2, plc)) // checks the geometry of segments
					{
						plc.mark(segIter1, sewedMark);
						plc.mark(segIter2, sewedMark);
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
		if (plc.is_sewable<2>(dIter->first, dIter->second))
		{
			plc.sew<2>(dIter->first, dIter->second);
			k++;	
		}

	cout << "\nNumber of sewable facets: " << k;
	cout << "\nNumber of pairs to be sewed: " << twoCellsToBeSewed.size() << "\n";
}



void writePLYOutput(LCC &lcc, string fileName)
{
	p_ply lccOutputPLY;

	if ((lccOutputPLY = ply_create(filename.c_str(), PLY_ASCII, NULL, 0, NULL)) == NULL)
	{
		cout << "\nCannot open file for writing!!";
		exit(0);
	}

	// count number of vertices and faces in LCC
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
			ply_write(lccOutputPLY, lcc.info<0>(pointInFaceIter)); 
	}

	ply_close(lccOutputPLY);			
}


// computes delaunay tetrahedralization
void computeDelaunayTetrahedralization()
{	
	vector <pair <CGALPoint, DartHandle> > lccVertexVector;

	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
		lccVertexVector.push_back(plc.point(pIter), pIter);

	DT.insert(lccVertexVector.begin(), lccVertexVector.end());

	cout << "\nDelaunay tetrahedralization computed!!";
	cout << "\nNumber of vertices in Delaunay tetrahedralization:" << DT.number_of_vertices();
	cout << "\nNumber of tetrahedrons in Delaunay tetrahedralization:" << DT.number_of_cells();
}



///////////////////////////////////////////// Segment recovery //////////////////////////////////////////////////////////////


void formMissingSegmentsQueue(vector<DartHandle> &missingSegmentQueue)
{
	missingSegmentQueue.clear();

	for (LCC::One_dart_per_cell_range<1>::iterator segmentIter = plc.one_dart_per_cell<1>().begin(), segmentIterEnd = plc.one_dart_per_cell<1>().end(); segmentIter != segmentIterEnd; segmentIter++)
	{
		CGALPoint p1 = plc.point(segmentIter);
		CGALPoint p2 = plc.point(plc.beta(segmentIter, 1));

		Delaunay::Cell_handle c;
		int i, j;

		if (DT.is_vertex(p1, vh1))
			if (DT.is_vertex(p2, vh2))
				if (!DT.is_edge(vh1, vh2, c, i, j))
					missingSegmentQueue.push_back(segmentIter);
	}

	cout << "\nTotal number of missing constraint segments:" << missingSegmentQueue.size() << "\n";


	return;
}

unsigned int computeCircumradius(Point &A, Point &B, Point &encroachingCandidate)
{
	// computing circumradius of a triangle:
	
	float a, b, c;
	float circumradius;
	a = sqrt(pow((A.x() - encroachingCandidate.x()), 2)+ pow((A.y() - encroachingCandidate.y()), 2) + pow((A.z() - encroachingCandidate.z()), 2));
	b = sqrt(pow((B.x() - encroachingCandidate.x()), 2)+ pow((B.y() - encroachingCandidate.y()), 2) + pow((B.z() - encroachingCandidate.z()), 2));
	c = sqrt(pow((A.x() - B.x()), 2)+ pow((A.y() - B.y()), 2) + pow((A.z() - B.z()), 2));

	return circumradius = (a * b * c) / sqrt((a + b + c) * (b + c - a) * (c + a - b) * (a + b - c));
}

void computeReferencePoint(Point *refPoint, DartHandle missingSegmentHandle)
{

	Point &A = plc.point(missingSegmentHandle);
	Point &B = plc.point(plc.beta(missingSegmentHandle, 1));

	float missingSegmentLength = sqrt(pow((A.x() - B.x()), 2) + pow((A.y() - B.y()), 2) + pow((A.z() - B.z()), 2));
	float sphereRadius = (missingSegmentLength / 2.0);
	CGALPoint sphereCenter = Point((A.x() + B.x()) / 2.0,  (A.y() + B.y()) / 2.0, (A.z() + B.z()) / 2.0);

	Sphere smallestCircumsphere(sphereCenter, pow(sphereRadius, 2));
	
	float encroachingCandidateDistance;
	vector<float> circumradiusMap; // value is computed only for those points which are encroaching
	
	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>.end(); pIter != pIterEnd; pIter++)
	{
		CGALPoint candidatePoint = plc.point(pIter);
		
		if (candidatePoint != A && candidatePoint != B)
		{
			encroachingCandidateDistance = sqrt(pow(sphereCenter.x() - candidatePoint.x(), 2) + pow(sphereCenter.y() - candidatePoint.y(), 2) + pow(sphereCenter.z() - candidatePoint.z(), 2));
	
			if (encroachingCandidateDistance <= sqrt(smallestCircumsphere.squared_radius())) // encroaching vertex
				circumradiusMap.push_back(computeCircumradius(A, B, candidatePoint));
			else // not encroaching
				circumradiusMap.push_back(INVALID_VALUE);
		
		}
		else
			continue; // skip
	}

	// select the one with maximum circumradius.
	float maxCircumradius = INVALID_VALUE;
	size_t vertexId = 0;

	for (LCC::One_dart_per_cell_range<0>::iterator pointIter = plc.one_dart_per_cell<0>().begin(), pointIterEnd = plc.one_dart_per_cell<0>.end(); pointIter != pointIterEnd; pointIter++, vertexId++)
		if (circumradiusMap[vertexId] != INVALID_VALUE)
			if (maxCircumradius < circumradiusMap[vertexId])
			{
				maxCircumradius = circumradiusMap[vertexId];
				refPoint = &(plc.point(pointIter));
			}	

	return;
}

float dotProduct (DartHandle segment1Handle, DartHandle segment2Handle)
{
	CGALPoint segment1Vertex[2];
	CGALPoint segment2Vertex[2];

	segment1Vertex[0] = plc.point(segment1Handle);
	segment1Vertex[1] = plc.point(plc.beta(segment1Handle, 1));

	segment2Vertex[0] = plc.point(segment2Handle);
	segment2Vertex[1] = plc.point(plc.beta(segment2Handle, 1));

	CGALPoint vector1 = CGALPoint(segment1Vertex[0].x() - segment1Vertex[1].x(), segment1Vertex[0].y() - segment1Vertex[1].y(), segment1Vertex[0].z() - segment1Vertex[1].z());
	CGALPoint vector2 = CGALPoint(segment2Vertex[0].x() - segment2Vertex[1].x(), segment2Vertex[0].y() - segment2Vertex[1].y(), segment2Vertex[0].z() - segment2Vertex[1].z());


	float v1Dotv2 = vector1.x() * vector2.x() + vector1.y() * vector2.y() + vector1.z() * vector2.z();

	return v1Dotv2;

}

float vectorMagnitude(DartHandle inputSegmentHandle)
{
	CGALPoint segmentVertices[2];

	segmentVertices[0] = plc.point(inputSegmentHandle);
	segmentVertices[1] = plc.point(plc.beta(inputSegmentHandle, 1));	

	CGALPoint vector = Point(segmentVertices[0].x() - segmentVertices[1].x(), segmentVertices[0].y() - segmentVertices[1].y(), segmentVertices[0].z() - segmentVertices[1].z());

	float vectorMagnitude = sqrtf(powf(vector.x(), 2.0) + powf(vector.y(), 2.0) + powf(vector.z(), 2.0));

	return vectorMagnitude;
}


float computeAngleBetweenSegments(DartHandle segment1Handle, DartHandle segment2Handle)
{

	float angle = acosf(dotProduct(segment1Handle, segment2Handle) / (vectorMagnitude(segment1Handle) * vectorMagnitude(segment2Handle)));

	return (angle * 180.0f / Pi); // conversion to degrees
}

bool isVertexAcute(DartHandle inputPointHandle)
{
	// Determine segment-pair(involving A) 
	vector<DartHandle> incidentOnInputPoint;

	 
	for (LCC::One_dart_per_incident_cell_range<1, 0>::iterator incidentSegmentIter = plc.one_dart_per_incident_cell<1, 0>(inputPointHandle).begin(), incidentSegmentIterEnd = plc.one_dart_per_incident_cell<1, 0>(inputPointHandle).end(); incidentSegmentIter++)
		incidentOnInputPoint.push_back(incidentSegmentIter);

	// Compute angle between all possible pairs(NAIVE SOLUTION)
	for (vector<DartHandle>::iterator segIter1 = incidentOnInputPoint.begin(); segIter1 != incidentOnInputPoint.end(); segIter1++) 
		for (vector<DartHandle>::iterator segIter2 = incidentOnInputPoint.begin(); *segIter1 != *segIter2 && segIter2 != incidentOnInputPoint.end(); segIter2++)
			if (computeAngleBetweenSegments(segIter1, segIter2) < 90.0f)
				return true;

	return false; 
}


unsigned int determineSegmentType(DartHandle missingSegmentHandle)
{

	Point &A = plc.point(missingSegmentHandle);
	Point &B = plc.point(plc.beta(missingSegmentHandle, 1));

	bool vertexAIsAcute = isVertexAcute(missingSegmentHandle);
	bool vertexBIsAcute = isVertexAcute(plc.beta(missingSegmentHandle, 1));

	if (!vertexAIsAcute && !vertexBIsAcute)	
		return 1;
	
	else if (vertexAIsAcute != vertexBIsAcute) // effectively XOR
		return 2;
	else
		return 3;
}

float computeSegmentLength(Point &A, Point &B)
{
	float sLength;
	sLength = sqrt(pow(A.x() - B.x(), 2) + pow(A.y() - B.y(), 2) + pow(A.z() - B.z(), 2));
	return sLength;
}

void updatePLCAndDT(Point &v, DartHandle missingSegmentHandle)
{

	// update PLC
	insert_point_in_cell<1>(missingSegmentHandle, v);
	// update DT
	computeDelaunayTetrahedralization(); 
}


void splitMissingSegment(DartHandle missingSegmentHandle)
{

	CGALPoint vb, refPoint;
	CGALPoint sphereCenter;
	float sphereRadius;
	unsigned int segmentType;
	segmentType = determineSegmentType(missingSegmentHandle);

	computeReferencePoint(&refPoint, missingSegmentId);
	
	CGALPoint A = plc.point(missingSegmentHandle);
	CGALPoint B = plc.point(plc.beta(missingSegmentHandle, 1));
	
	float AP, PB, AB;

	CGALPoint v;

	if (segmentType == 1)
	{
		AP = computeSegmentLength(A, refPoint);
		AB = computeSegmentLength(A, B);
		
		if (AP < 0.5f * AB)
		{
			sphereCenter = A;
			sphereRadius = AP;
		}

		else if (PB < 0.5 * AB)
		{
			sphereCenter = B;
			sphereRadius = PB;
		}

		else
		{
			sphereCenter = A;
			sphereRadius = 0.5 * AB; 
		}	

		CGALSphericalPoint sphericalSphereCenter = CGALSphericalPoint(sphereCenter.x(), sphereCenter.y(), sphereCenter.z());
		CGALSphericalSphere s = CGALSphericalSphere(sphericalSphereCenter, pow(sphereRadius, 2));

		CGALSphericalPoint p1 = CGALSphericalPoint(A.x(), A.y(), A.z());
		CGALSphericalPoint p2 = CGALSphericalPoint(B.x(), B.y(), B.z());

		CGALSegment seg(p1, p2);
		CGALSphericalSegment lineArc(seg);
		vector<Object> intersections;
		
		CGAL::intersection(lineArc, s, back_inserter(intersections));
		
		if (intersections.size() > 0)
		{	v = object_cast<CGALPoint>(intersections.back());
			intersections.pop_back();
		}

		else
		{
			cout << "Case 1: Sphere and segment do not intersect";
			exit(0);
		}
	}

	else if (segmentType == 2)
	{
		DartHandle acuteParentHandle; 
		DartHandle ApB[1]; 
		float vbLength;
		DartHandle AHandle = missingSegmentHandle;
		DartHandle BHandle = plc.beta(missingSegmentHandle, 1);

		if (isVertexAcute(AHandle) && !isVertexAcute(BHandle))
			acuteParentHandle = AHandle;
	        else if (isVertexAcute(BHandle) && !isVertexAcute(AHandle))
		{
			acuteParentHandle = BHandle;
			// swap A <-> B 
			unsigned int t;
			t = A;
			A = B;
			B = t;
		}
		
		
		// Segment calculations
		
		ApB[0] = acuteParentHandle;
		ApB[1] = BHandle;
		
	
		CGALSphericalPoint p1(plc.point(acuteParentHandle).x(), plc.point(acuteParentHandle).y(), plc.point(acuteParentHandle).z());
		CGALSphericalPoint p2(plc.Point(BHandle).x(), plc.point(BHandle).y(), plc.point(BHandle).z());
		
		CGALSegment seg(p1, p2);
		CGALSphericalSegment lineArc(seg);
		
		/// Sphere calculations
		CGALSphericalPoint acuteParent(plc.point(acuteParentHandle).x(), plc.point(acuteParentHandle).y(), plc.point(acuteParentHandle).z());
		
		CGALPoint acuteParentLinearField(plc.point(acuteParentHandle).x(), plc.point(acuteParentHandle).y(), plc.point(acuteParentHandle).z());
	
		float ApRefPointLength = computeSegmentLength(acuteParentLinearField, refPoint);

		CGALSphericalSphere s(acuteParent, pow(ApRefPointLength, 2));

		vector<Object> intersections;

		CGAL::intersection(lineArc, s, back_inserter(intersections));
	
		if (intersections.size() > 0)
		{	
			v = object_cast<Point>(intersections.back());
			intersections.pop_back();
		}

		else
		{
			cout << "\nCase 2: Sphere and segment do not intersect";
			exit(0);
		}
		
		float vrefpointLength = computeSegmentLength(v, refPoint);
		
		float acuteparentALength;
			
		if (vbLength < vrefpointLength) // v was rejected
		{
			SphericalPoint sphereCenter(acuteParent);
			unsigned int avLength = computeSegmentLength(plcVertices[A].first, v);
			if (vrefpointLength < 0.5 * avLength)
				{
					Point temp(acuteParentLinearField.x(), acuteParentLinearField.y(), acuteParentLinearField.z());
					acuteparentALength = computeSegmentLength(temp, plcVertices[A].first);
					sphereRadius = acuteparentALength + avLength - vrefpointLength;	
				}
			else
				sphereRadius = acuteparentALength + 0.5 * avLength;
		        
			s = SphericalSphere(sphereCenter, pow(sphereRadius, 2));
		
			CGAL::intersection(lineArc, s, back_inserter(intersections));

			if (intersections.size() > 0)
			{
				v = object_cast<Point>(intersections.back());
				intersections.pop_back();
			}
			else
			{
				cout << "Case 3: Sphere and segment do not intersect";
				exit(0);
			}
		}
	}

	else if (segmentType == 3) // meaning both A & B are acute
	{
		CGALPoint newPoint;
		
		float x = (A.x() + B.x()) / 2.0;
		float y = (A.y() + B.y()) / 2.0;
		float z = (A.z() + B.z()) / 2.0;
		newPoint = CGALPoint(x, y, z); 
		
		v = newPoint;
	}

	// update plc and DT
	updatePLCAndDT(v, missingSegmentHandle);
	
	return;
}


void recoverConstraintSegments()
{
	// I/P: plc1, DT1
	// O/P: plc2, DT2
	
	vector<DartHandle> missingSegmentQueue;
	DartHandle missingSegment;

	formMissingSegmentsQueue(missingSegmentQueue);
	int i = 0;
	while (missingSegmentQueue.size() != 0)
	{
		missingSegment = missingSegmentQueue.back();
		missingSegmentQueue.pop_back();
		splitMissingSegment(missingSegment);
		formMissingSegmentsQueue(missingSegmentQueue);
		i++;
	}
	
	cout << "\nIn the loop " << i << " number of times" << "\n";
	
	return;
}

////////////////////////////////////////////// Local Degeneracy Removal begin ///////////////////////////////////////////////////

class DegenerateVertexSetCandidate
{
	public:
		DartHandle pointHandles[5];
};


// returns true if  given vertices are co-spherical
bool areCospherical(DegenerateVertexSetCandidate degenSet)
{	
	CGALPoint p[5];
	
	for (unsigned int i = 0; i < 5; i++)	
		p[i] = plc.point(degenSet[i]);

	if (side_of_bounded_sphere(p[0], p[1], p[2], p[3], p[4]) == ON_BOUNDARY)
		return true;
	else
		return false;
}

void addLocalDegeneraciesToQueue(vector<DegenerateVertexSetCandidate> &localDegeneracySet)
{	
	DegenerateVertexSetCandidate degenerateSetCandidate;
	
	for (Delaunay::Finite_cells_iterator cellIter = DT.finite_cells_begin(); cellIter != DT.finite_cells_end(); cellIter++)
	{
		for (unsigned int n = 0; n < 3; n++)
			degenerateSetCandidate.pointIds[n] = (cellIter->vertex(n))->info(); // info structure contains pointIds

		for (unsigned int j = 0; j < 4; j++)
			{
				VertexHandle vh = DT.mirror_vertex(cellIter, j);	
				degenerateSetCandidate.pointIds[4] = vh->info();		
						
				if (areCospherical(degenerateSetCandidate))
				{
					localDegeneracySet.push_back(degenerateSetCandidate);	
				}
			}
	}
}

			
// perturbation should be such that it does not make PLC inconsistent
bool isVertexPerturbable(unsigned pointId)
{
	bool perturbable = false;
	
		// if the vertex is collinear
		// if the vertex is coplanar
			// Yes
		// else
			// No 	 		

	// for all edges incident on vertex 'pointId'
		// check is eny 2 edges are collinear
		// if yes, then check coplanarity
		// else, return false
		// check coplanarity:
		// for all faces containing vertex 'pointId'
		// check if there are atleast a pair of faces which are coplanar
		// If yes, return true
		// else, return false


	return perturbable;
}

bool isVertexSegmentSafePertubable(unsigned int pointId)
{
	bool segmentSafePerturbable = false;

	return segmentSafePerturbable;
}

bool isDegeneracyRemovable(DegenerateVertexSetCandidate degenCandidate)
{
	bool removable = false;
	


	return removable;
}

void perturbRemove(unsigned int pointId, vector<DegenerateVertexSetCandidate>& localDegeneracySet)
{
	// DO NOTHING FOR NOW
}


bool areAffinelyIndependent(DegenerateVertexSetCandidate degenSet)
{

	CGALPoint v1 = plc.point(degenSet[pointId].pointHandles[0]);
	CGALPoint v2 = plc.point(degenSet[pointId].pointHandles[1]);
	CGALPoint v3 = plc.point(degenSet[pointId].pointHandles[2]);
	CGALPoint v4 = plc.point(degenSet[pointId].pointHandles[3]);
	CGALPoint v5 = plc.point(degenSet[pointId].pointHandles[4]);

	if (coplanar(v1, v2, v3, v4) || coplanar(v2, v3, v4, v5) || coplanar(v3, v4, v5, v1) || coplanar(v4, v5, v1, v2) || coplanar(v5, v1, v2, v3))
		return false;

	return true; 
}

void computeBreakPoint(CGALPoint &vb, unsigned int pointId, vector<DegenerateVertexSetCandidate> &localDegeneracySet)
{

	CGALPoint v1 = plc.point(localDegeneracySet[pointId].pointHandles[0]);
	CGALPoint v2 = plc.point(localDegeneracySet[pointId].pointHandles[1]);
	CGALPoint v3 = plc.point(localDegeneracySet[pointId].pointHandles[2]);
	CGALPoint v4 = plc.point(localDegeneracySet[pointId].pointHandles[3]);
	CGALPoint v5 = plc.point(localDegeneracySet[pointId].pointHandles[4]);

	if (areAffinelyIndependent(localDegeneracySet[pointId]))
	{
		CGALSphere s(v1, v2, v3, v4);
		vb = s.center();
	}	
	else
	{
		CGALCircle c;

		if (coplanar(v1, v2, v3, v4))	
			c = CGALCircle(v1, v2, v3);

		if (coplanar(v2, v3, v4, v5))	
			c = CGALCircle(v2, v3, v4);

		if (coplanar(v3, v4, v5, v1))	
			c = CGALCircle(v3, v4, v5);

		if (coplanar(v4, v5, v1, v2))	
			c = CGALCircle(v4, v5, v1);

		if (coplanar(v5, v1, v2, v3))	
			c = CGALCircle(v5, v1, v2);

		vb = c.center(); 
	}
}

bool isVertexEncroachingSegment(CGALPoint vb, DartHandle segmentHandle)
{

	CGALPoint endpoint1(plc.point(segmentHandle));
	CGALPoint endpoint2(plc.point(plc.beta(segmentHandle, 1)));

	float x = (endpoint1.x() + endpoint2.x()) / 2.0;
	float y = (endpoint1.y() + endpoint2.y()) / 2.0;
	float z = (endpoint1.z() + endpoint2.z()) / 2.0;

	CGALPoint sphereCenter(x, y, z);

	float sphereRadius = sqrt(pow(sphereCenter.x() - endpoint1.x(), 2) + pow(sphereCenter.y() - endpoint1.y(), 2) + pow(sphereCenter.z() - endpoint1.z(), 2));

	float centerToBreakingPointDistance = sqrt(pow(sphereCenter.x() - vb.x(), 2) + pow(sphereCenter.y() - vb.y(), 2) + pow(sphereCenter.z() - vb.z(), 2));


	if (sphereRadius <= centerToBreakingPointDistance)
		return true;
	else
		return false;
		
}

bool isVertexEncroachingFace(Point vb, DartHandle facetHandle)
{
	
	CGALPoint v1(plc.point(facetHandle));
	CGALPoint v2(plc.point(plc.beta(facetHandle, 1)));
	CGALPoint v3(plc.point(plc.beta(facetHandle, 1, 1)));

	CGALCircle c(v1, v2, v3);
	CGALPoint circleCenter(c.center());
	float centerToBreakingPointDistance = sqrt(pow(circleCenter.x() - vb.x(), 2) + pow(circleCenter.y() - vb.y(), 2) + pow(circleCenter.z() - vb.z(), 2));

	if (sqrtf(c.squared_radius()) < centerToBreakingPointDistance)
		return true;

	return false;
}


bool isEncroachingPLC(CGALPoint vb)
{
	// encroaches any segment
	for (LCC::One_dart_per_cell_range<1>::iterator segmentIter = plc.one_dart_per_cell<1>().begin(), segmentIterEnd = plc.one_dart_per_cell<1>().end(); segIter != segIterEnd; segIter++)
		if (isVertexEncroachingSegment(vb, segIter))
			return true;

	// encroaches any face
	for (unsigned int k = 0; k < plcFaces.size(); k++)
		if (isVertexEncroachingFace(vb, k))
			return true;

	return false;
}

void boundaryProtection(CGALPoint vb)
{

	for (LCC::One_dart_per_cell_range<1>::iterator segmentIter = plc.one_dart_per_cell<1>().begin(), segmentIterEnd = plc.one_dart_per_cell<1>().end(); segIter != segIterEnd; segIter++)
		if (isVertexEncroachingSegment(vb, segIter))
		{
			// circumcenter of segment 'n':
			CGALPoint endpoint1(segmentIter);
		        CGALPoint endpoint2(plc.beta(segmentIter, 1));	
			float x = (endpoint1.x() + endpoint2.x()) / 2.0;
			float y = (endpoint1.y() + endpoint2.y()) / 2.0;
			float z = (endpoint1.z() + endpoint2.z()) / 2.0;

			Point sphereCenter(x, y, z);

			updatePLCAndDT(sphereCenter, segmentIter);
		}
	}
	
}


// removes local degeneracies from Delaunay tetrahedralization
void removeLocalDegeneracies()
{

// I/P: plcVertices1, plcFaces1, plcSegments1, DT1
// O/P: plcVertices2, plcFaces2, plcSegments2, DT2

	cout << "\nStarting local degeneracy removal...";

	// compute all local degeneracies in DT and add them to Q
	vector<DegenerateVertexSetCandidate> localDegeneracySet;
	addLocalDegeneraciesToQueue(localDegeneracySet);

	// repeat while Q != NULL	
	while (localDegeneracySet.size() != 0)
	{	
		for (unsigned int n = 0; n < localDegeneracySet.size(); n++)
		{
			if (isDegeneracyRemovable(localDegeneracySet[n]))
				perturbRemove(n, localDegeneracySet); 
			else
			{
				CGALPoint vb;	
				computeBreakPoint(vb, n, localDegeneracySet);
				if (isEncroachingPLC(vb))
					boundaryProtection(vb);
				else
					plc.add_vertex(vb);
			}						
		}
	}
	cout << "\nLocal degeneracy removal completed";
}

///////////////////////////////////////////////////Local Degeneracy Removal Ends//////////////////////////////////////////////////////



/////////////////////////////////////////////// Facet recovery starts ///////////////////////////////////////////////////////////

void formMissingSubfaceQueue(vector<DartHandle> &missingSubfacesQueue)
{

	Delaunay::Cell_handle ch;
	int i, j, k;
	
	Delaunay::Vertex v1, v2, v3;
	VertexHandle vh1, vh2, vh3;

	for (LCC::One_dart_per_cell_range<2>::iterator facetIter = plc.one_dart_per_cell<2>().begin(), facetIterEnd = plc.one_dart_per_cell<2>().end(); facetIter != facetIterEnd; facetIter++)
	{
		DT.is_vertex(plc.point(facetIter), vh1);
		DT.is_vertex(plc.point(plc.beta(facetIter, 1)), vh2);
		DT.is_vertex(plc.point(plc.beta(facetIter, 1, 1)), vh3);

		if (DT.is_facet(vh1, vh2, vh3, ch, i, j, k))
			continue;
		else
			missingSubfacesQueue.push_back(facetIter);
	}
}


void copyInfoFromDTToLCC(Delaunay dt, LCC& linearCellComplex, map <CellHandle, DartHandle> *dtCellToLCCCellMap)
{

	for (Delaunay::Finite_cells_iterator delaunayCellIter = dt.finite_cells_begin(); delaunayCellIter != dt.finite_cells_end(); delaunayCellIter++)
	{
		DartHandle lccCellHandle = (*dtCellToLCCCellMap)[delaunayCellIter];

		for (unsigned int n = 0; n < 4; n++)
		{
			CGALPoint delaunayPoint = (delaunayCellIter->vertex(n))->point();

			for (LCC::One_dart_per_incident_cell_range<0, 3>::iterator lccPointIter = linearCellComplex.one_dart_per_incident_cell<0, 3>(lccCellHandle).begin(); lccPointIter != linearCellComplex.one_dart_per_incident_cell<0, 3>(lccCellHandle).end(); lccPointIter++)
			{
				CGALPoint lccPoint = linearCellComplex.point(lccPointIter);
				if (delaunayPoint == lccPoint)			
				{
					linearCellComplex.info<0>(lccPointIter) = (delaunayCellIter->vertex(n))->info();
					break;
				}
			}
		}
	}

}


void createEquivalentTetrahedralization()
{
	
	map<CellHandle, DartHandle> *dtVolumeToLccDartMap;
	import_from_triangulation_3(cdtMesh, DT, dtVolumeToLccDartMap);

	copyInfoFromDTToLCC(DT, cdtMesh, dtVolumeToLccDartMap);
}


void formCavity(vector<DartHandle> *cavity, DartHandle missingSubfaceHandle, vector<DartHandle>& lcc3CellsToBeRemoved)
{
	vector<DartHandle> intersectingTets;				
	Point pTet[4], pTri[3];
	
	pTri[0] = plc.point(missingSubfaceHandle);
	pTri[1] = plc.point(plc.beta(missingSubfaceHandle, 1));
	pTri[2] = plc.point(plc.beta(missingSubfaceHandle, 1, 1));
	 
	// determine the intersecting 3-cells
	unsigned int i = 0;
	for (LCC::One_dart_per_cell_range<3>::iterator cellIter = cdtMesh.one_dart_per_cell<3>().begin(), cellIterEnd = cdtMesh.one_dart_per_cell<3>().end(); cellIter != cellIterEnd; cellIter++)
	{
			
		i = 0;

		for (LCC::One_dart_per_incident_cell_range<0, 3>::iterator vertexIter = cdtMesh.one_dart_per_incident_cell<0, 3>(cellIter).begin(), vertexIterEnd = cdtMesh.one_dart_per_incident_cell<0, 3>(cellIter).end(); vertexIter != vertexIterEnd; vertexIter++)
			pTet[i++] = cdtMesh.point(vertexIter);
		
		
		CGALTetrahedron CGALTet(pTet[0], pTet[1], pTet[2], pTet[3]);
		
		CGALTriangle CGALTri(pTri[0], pTri[1], pTri[2]);
	
		if (do_intersect(CGALTri, CGALTet))
		{
			intersectingTets.push_back(cellIter);
			lcc3CellsToBeRemoved.push_back(cellIter);
		}
		else
			continue;
	}

	// Global cavity formation
	map<DartHandle, unsigned int> facetVisitCounterMap;

 	// Initialize facetVisitCounter
	for (LCC::One_dart_per_cell_range<2>::iterator faceIter = cdtMesh.one_dart_per_cell<2>().begin(), faceIterEnd = cdtMesh.one_dart_per_cell<2>().end(); facetIter != facetIterEnd; faceIter++)
		facetVisitCounterMap.insert(pair<DartHandle, unsigned int>(faceIter, 2)); // TODO: will it have 1 dart for a facet or it will contain pair for darts(in opposite direction for each face) ???

	DartHandle tempTetId;
	// Compute facetVisitCounter
	for (unsigned int n = 0; n < intersectingTets.size(); n++)
	{
		tempTetId  = intersectingTets.back();
		intersectingTets.pop_back();
		
		for (LCC::One_dart_per_incident_cell_range<2 ,3>::iterator fIter = cdtMesh.one_dart_per_incident_cell<2, 3>(tempTetId).begin(), fIterEnd = cdtMesh.one_dart_per_incident_cell<2, 3>(tempTetId).end(); fIter != fIterEnd; fIter++)
			facetVisitCounterMap[fIter] = facetVisitCounterMap[fIter] - 1; // TODO: Doubtful
										   
	}	

	// Determine globalCavity using faceVisitCounter	
	vector<DartHandle> globalCavity;	
	for (map<DartHandle, unsigned int>::iterator iter = facetVisitCounterMap.begin(), iterEnd = facetVisitCounterMap.end(); iter != iterEnd; iter++)
		if (facetVisitCounterMap[iter->first] == 1)
			globalCavity.push_back(iter->first);
		else
			continue;	
		

	// Partition global cavity to upper and lower cavity
		// For each face in global cavity determine position of a point on its surface wrt. missingSubface
		// If position of this point is above this face put this face in upper cavity otherwise in bottom cavity
	
	CGALPoint v1 = plc.point(missingSegmentHandle); 
	CGALPoint v2 = plc.point(plc.beta(missingSegmentHandle, 1));
	CGALPoint v3 = plc.point(plc.beta(missingSegmentHandle, 1, 1));	

	// generating random points
	Random rnd1(30.0), rnd2(10.0);

	float s = rnd1.uniform_real<float>(); 
	float t = rnd2.uniform_real<float>();
	CGALPoint points[3];
	
	float xRand, yRand, zRand;


	for (unsigned int g = 0; g < globalCavity.size(); g++)
	{
		n = 0;		
		
		for (LCC::One_dart_per_incident_cell_range<0, 2>::iterator vertexIter = cdtMesh.one_dart_per_incident_cell<0, 2>(globalCavity[g]).begin(); vertexIter != cdtMesh.one_dart_per_incident_cell<0, 2>(globalCavity[g]).end(); vertexIter++)
			points[n++] = cdtMesh.point(vertexIter);
		
		xRand = (1.0 - s - t) * points[0].x() + s * points[1].x() + t * points[2].x(); // parametric representation of a point on the plane
		yRand = (1.0 - s - t) * points[0].y() + s * points[1].y() + t * points[2].y();	
		zRand = (1.0 - s - t) * points[0].z() + s * points[1].z() + t * points[2].z();
	
		CGALPoint randomPointOnTheFacet(xRand, yRand, zRand); // random point inside triangle 'g' of global cavity
			
		while (orientation(v1, v2, v3, randomPointOnTheFacet) == COPLANAR)
		{
			s = rnd1.uniform_real<float>();
			t = rnd2.uniform_real<float>();
			
			xRand = (1.0 - s - t) * points[0].x() + s * points[1].x() + t * points[2].x();
			yRand = (1.0 - s - t) * points[0].y() + s * points[1].y() + t * points[2].y();	
			zRand = (1.0 - s - t) * points[0].z() + s * points[1].z() + t * points[2].z();

			randomPointOnTheFacet = CGALPoint(xRand, yRand, zRand); // re-initialize
		}

		if (orientation(v1, v2, v3, randomPointOnTheFacet) == CGAL::POSITIVE)
			cavity[0].push_back(globalCavity[g]);
		else // case of coplanarity already removed
			cavity[1].push_back(globalCavity[g]);
		
	}	


}

bool isStronglyDelaunay(DartHandle facetHandle, vector<DartHandle> cavityVerticesSet, LCC cavityLCC)
{
	// tests whether facet/triangle pointed to by facetHandle is strongly Delaunay wrt. vertices in cavityVerticesSet

	// compute Delunay tetrahedralization of vertices
		// check if the face is there in that DT 
		// If yes then check if the enclosing sphere has any other vertex on its surface
			// If yes, then it is Delaunay but not strongly Delaunay, return false
			// If no, then it is strongly Delaunay, return true
		// If no, then it is not even Delaunay, return false 	
		
	Delaunay tempDT;
	vector<CGALPoint> cavityVertices;

	for (vector<DartHandle>::iterator vertexIter = cavityVerticesSet.begin(), vertexIter = cavityVerticesSet.end(); vertexIter != vertexIterEnd; vertexIter++)
		cavityVertices.push_back(cavityLCC.point(vertexIter));


	// Compute DT
	tempDT.insert(cavityVertices.begin(), cavityVertices.end());
	
	// test strong Delaunay criteria
	VertexHandle vh[3];
	int i, j, k, h = 0;
	CellHandle ch;

	for (LCC::One_dart_per_incident_cell_range<0, 2>::iterator iter = cavityLCC.one_dart_per_incident_cell<0, 2>(facetHandle).begin(), iterEnd != cdtMesh.one_dart_per_incident_cell<0, 2>(facetHandle).end(); iter != iterEnd; iter++)
		tempDT.is_vertex(cavityLCC.point(iter), vh[h++]);


	if (tempDT.is_facet(vh[0], vh[1], vh[2], ch, i, j, k))
	{
		CGALSphere s(vh[0]->point(), vh[1]->point(), vh[2]->point(), ((*ch).vertex(6 - i - j - k))->point());

		for (vector<Point>::iterator vIter = cavityVertices.begin(); vIter != cavityVertices.end(); vIter++)
		{	
			if (s.has_on_boundary(*vIter))
				return false;	
		}

		return true;
	}

	return false;

}



int locateFacetInCavity(DartHandle nonStronglyDelaunayFacet, vector<DartHandle> cavity)
{

	//if (facetPosition = locateFacetInCavity(tempNonStronglyDelaunayFacet, cavity[i])) // found, 
	int facetLocation = 0;
	for (vector<DartHandle>::iterator facetIter = cavity.begin(); facetIter != cavity.end(); facetIter++)
	{
		if (nonStronglyDelaunayFacet == *facetIter)
			return facetLocation;
		else
			facetLocation++;
	}
	return -1;
}


bool isCellOutsideCavity(DartHandle cellHandle, vector<DartHandle> cavity)
{


	// Approach 1:
		// Convert tet & cavity to Nef_polyhedra
		// Determine intersection
			// If intersection is tet itself
				// return false
			// else, 
				// return true	   
				
	// Approach 2:
		// Take a random ray out of barycenter of tetrahedron
		// Count number of intersections of the ray with polyhedron
		// If number of intersections:
			// Even:
				// return true
			// Odd:
				// return false 		
				
		Point p[4];

		unsigned int i = 0;
		for (lcc::One_dart_per_incident_cell_range<0, 3>::iterator vertexIter = cdtMesh.one_dart_per_incident_cell<0, 3>(cellHandle).begin(); vertexIter != cdtMesh.one_dart_per_incident_cell<0, 3>(cellHandle).end(); vertexIter++)
			p[i++] = plcVertices[cdtMesh.info<0>(vertexIter)].first;
		
		float x = 0.0, y = 0.0, z = 0.0;

		for (unsigned int n = 0; n < 4; n++)
		{
			x += p[n].x();
			y += p[n].y();
			z += p[n].z();
		}

			x /= 4.0;
			y /= 4.0;
			z /= 4.0;

		Point tetBarycenter(x, y, z);
	
		// generate a random ray
		
		Random rndGenerator1(50.0); // some numbers given as seed
		Random rndGenerator2(60.0);
		Random rndGenerator3(70.0);

		Point randomPoint(rndGenerator1.uniform_real<float>(), rndGenerator2.uniform_real<float>(), rndGenerator3.uniform_real<float>());

		Ray_3<K> randomRay(tetBarycenter, randomPoint);
		unsigned int intersectionCount = 0;

		for (vector<DartHandle>::iterator facetIter = cavity.begin(); facetIter != cavity.end(); facetIter++)
		{
			// compute its intersection with all faces of cavity
			// count number of intersections
			Point facetVertices[3];
			
			unsigned int k = 0;	
			for (lcc::One_dart_per_incident_cell_range<0, 2>::iterator pointIter = cdtMesh.one_dart_per_incident_cell<0, 2>(*facetIter).begin(); pointIter != cdtMesh.one_dart_per_incident_cell<0, 2>(*facetIter).end(); pointIter++)
				facetVertices[k++] = plcVertices[cdtMesh.info<0>(pointIter)].first;
		

			Triangle_3<K> cavityFacet(facetVertices[0], facetVertices[1], facetVertices[2]); 
		
			if (do_intersect(randomRay, cavityFacet))
				intersectionCount++;
			else
				continue;
		}
		
		if (intersectionCount % 2 == 0)
			return true;
		else
			return false;
}


bool areSameFacets(LCC lcc1, DartHandle d1, LCC lcc2, DartHandle d2)
{
	DartHandle d3 = lcc2.make_triangle(lcc1.point(d1), lcc1.point(lcc1.beta(d1, 1)), lcc1.point(lcc1.beta(d1, 1, 1))); // TODO: Fix this roundabout way of comparing facets from 2 LCCs
	bool areSameFacets = lcc2.are_facet_same_geometry(d3, d2);
	lcc2.remove_cell<2>(d3);	

	return areSameFacets;
}

DartHandle locateFacetInCavity(DartHandle facet, LCC cavityLCC)
{
	for (LCC::One_dart_per_cell<2>::iterator facetIter = cavityLCC.one_dart_per_cell<2>().begin(), facetIterEnd = cavityLCC.one_dart_per_cell<2>().end(); facetIter != facetIterEnd; facetIter++)
	{
		if (areSameFacets(facet, cdtMesh, facetIter, cavityLCC))
			return facetIter;
		else
			continue;
	}
}


void addFaceToLCC(CGALPoint p1, CGALPoint p2, CGALPoint p3, LCC &cavityLCC)
{
	DartHandle newFace = cavityLCC.make_triangle(p1, p2, p3);
	vector <pair <DartHandle, DartHandle> > twoCellsToBeSewed;
	int sewedMark = cavityLCC.get_new_mark();

	if (sewedMark == -1)
	{
		cout << "\nNo free mark available!!";
		exit(0);
	}

	// sew the triangle with existing faces in cavityLCC
	for (LCC::One_dart_per_incident_cell_range<1, 2>::iterator segIter1 = cavityLCC.one_dart_per_incident_cell(newFace).begin(), segIterEnd1 = cavityLCC.one_dart_per_incident_cell(newFace).end(); segIter1 != segIterEnd1; segIter1++)
	{
		if (!cavityLCC.is_marked(segIter1, sewedMark)) // not sewed till now
		{
			for (LCC::Dart_range::iterator segIter2 = cavityLCC.darts().begin(), segIterEnd2 = cavityLCC.darts().end(); segIter2 != segIterEnd2; segIter2++)
			{
				if (!cavityLCC.is_marked(segIter2, sewedMark) && cavityLCC.is_sewable<2>(segIter1, segIter2))
				{
					if (areGeometricallySameSegments(segIter1, segIter2, cavityLCC)) 
					{
						cavityLCC.mark(segIter1, sewedMark);
						cavityLCC.mark(segIter2, sewedMark);
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
		if (cavityLCC.is_sewable<2>(dIter->first, dIter->second))
		{
			cavityLCC.sew<2>(dIter->first, dIter->second);
			k++;	
		}

}



void cavityRetetrahedralization(vector <DartHandle>& cavity, vector<DartHandle>& lcc3CellsToBeRemoved)
{
	// cavity verification/expansion
		// for all faces 'f' in upper(or lower) cavity, 
			// Form a queue Q of all non-strongly Delaunay faces in cavity C
			// If f is not strongly Deluanay, 
				// Remove f from cavity list
				// For all tetrahedra t sharing f
					// for all other faces F of t F != f:
						// If F is not already in cavity 
							// add it to cavity  
						// Else 
							// remove it from cavity
						// Modify C and set of vertices V of cavity
			// repeat untill all faces of cavity are strongly Delaunay


	LCC cavityLCC;
	createLCCRepresentationForCavity(cavity, cavityLCC, cavityLCCTocdtMeshDartMap);



	vector<DartHandle> nonStronglyDelaunayFacesInCavity;
	vector<DartHandle> cavityVerticesSet; // V

	// compute set of non-strongly Delunay faces in cavity
	do
	{	
		for (LCC::One_dart_per_cell_range<0>::iterator vertexIter = cavityLCC.one_dart_per_cell<0>().begin(), vertexIterEnd = cavityLCC.one_dart_per_cell<0>().end(); vertextIter != vertexIterEnd; vertexIter++)
			cavityVerticesSet.push_back(vertexIter);

		for (LCC::One_dart_per_cell_range<2>::iterator facetIter = cavityLCC.one_dart_per_cell<2>().begin(), facetIterEnd = cavityLCC.one_dart_per_cell<2>().end(); facetIter != facetIterEnd; facetIter++)
		{
			if (isStronglyDelaunay(facetIter, cavityVerticesSet, cavityLCC))
		        	continue;
			else
				nonStronglyDelaunayFacesInCavity.push_back(facetIter);
		}
	
		// For each non strongly Delunay face
		// remove face from cavity[i]
		// find tets from cdtMesh sharing it:
		// for each such tet t:
			// for all facets F of t, F != f
			// if F is in cavity[i]
				// remove F from cavity[i]
			// else,
				// add F to cavity[i]
				
		for (unsigned int h = 0; h < nonStronglyDelaunayFacesInCavity.size(); h++)
		{
			DartHandle tempNonStronglyDelaunayFacet = nonStronglyDelaunayFacesInCavity.back();
			nonStronglyDelaunayFacesInCavity.pop_back();
		
			// find tet from cdtMesh sharing it:
			DartHandle cellToBeRemovedHandle = cavityLCCTocdtMeshDartMap[tempNonStronglyDelaunayFacet];
			lcc3CellsToBeRemoved.push_back(cellToBeRemovedHandle);

			for (LCC::One_dart_per_incident_cell_range<2, 3>::iterator faceIter = cdtMesh.one_dart_per_incident_cell<2, 3>(cellToBeRemovedHandle).begin(), faceIterEnd != cdtMesh.one_dart_per_incident_cell<2, 3>(cellToBeRemovedHandle).end(); facetIter != facetIterEnd; faceIter++) // for all faces of the removed cell
			{
				if (!areSameFacet(faceIter, tempNonStronglyDelaunayFacet))
				{
					if ((facetHandle = locateFacetInCavity(faceIter, cavityLCC)) != -1)
						cavityLCC.remove_cell<2>(facetHandle);
					else // add facet to the cavity
						//cavityLCC.push_back(faceIter);
					{
						CGALPoint p1 = plc.point(faceIter);
						CGALPoint p2 = plc.point(plc.beta(faceIter, 1));
						CGALPoint p3 = plc.point(plc.beta(faceIter, 1, 1));
						
						addFaceToLCC(p1, p2, p3, cavityLCC); 
						
					}
				}
			}
			// remove facet from cavity
			cavityLCC.remove_cell<2>(tempNonStronglyDelaunayFacet); 
		}
	}while(nonStronglyDelaunayFacesInCavity.size() != 0);

	
	// Remove cells from cdtMesh(creates space for retetrahedralization)
	for (vector<DartHandle>::iterator cIter = lcc3CellsToBeRemoved.begin(); cIter != lcc3CellsToBeRemoved.end(); cIter++)
		cdtMesh.remove_cell<3>(*cIter);
	



	// Cavity retetrahedralization:
		// Compute DT of vertices of input cavity
		// import it to lcc
			
			// create an unordered set of vertex info 
			// insert all vertices of set to DT
	cavityVerticesSet.clear();

	for (vector<DartHandle>::iterator facetIter = cavity.begin(); facetIter != cavity.end(); facetIter++)
		for (lcc::One_dart_per_incident_cell_range<0, 2>::iterator vertexIter = cdtMesh.one_dart_per_incident_cell<0, 2>(*facetIter).begin(); vertexIter != cdtMesh.one_dart_per_incident_cell<0, 2>(*facetIter).end(); vertexIter++)
			cavityVerticesSet.insert(cdtMesh.info<0>(vertexIter));		
	

	// create a vector of vertices of cavity
	vector <pair<Point, unsigned int> > cavityVertices;
	for (unordered_set<unsigned int>::iterator pointIdIter = cavityVerticesSet.begin(); pointIdIter != cavityVerticesSet.end(); pointIdIter++)
		cavityVertices.push_back(plcVertices[*pointIdIter]);

	
	//create cavity vertices
	Delaunay cavityDT;
	lcc cavityLCC;

	cavityDT.insert(cavityVertices.begin(), cavityVertices.end());

	map<Delaunay::Cell_handle, lcc::Dart_handle> *dtCellToLCCCellMap;
	import_from_triangulation_3(cavityLCC, cavityDT, dtCellToLCCCellMap);
	copyInfoFromDTToLCC(cavityDT, cavityLCC, dtCellToLCCCellMap);

	// Sew retetrahedralized cavity back to original cdtMesh
	// Add cavityLCC to cdtMesh
	Point_3<K> tempPt[4];
	unsigned int x;

	for (lcc::One_dart_per_cell_range<3>::iterator cIter = cavityLCC.one_dart_per_cell<3>().begin(); cIter != cavityLCC.one_dart_per_cell<3>().end(); cIter++)
	{
		x = 0;

		for (lcc::One_dart_per_incident_cell_range<0, 3>::iterator vIter = cavityLCC.one_dart_per_incident_cell<0, 3>(cIter).begin(); vIter != cavityLCC.one_dart_per_incident_cell<0, 3>(cIter).end(); vIter++)
			tempPt[x++] = cdtMesh.point(vIter); 			
			
		cdtMesh.make_tetrahedron(tempPt[0], tempPt[1], tempPt[2], tempPt[3]);
	}

	cdtMesh.sew3_same_facets();
		
}


// recovers the constraint faces
void recoverConstraintFaces()
{
	 			
        vector<DartHandle> missingSubfacesQueue;
	formMissingSubfaceQueue(missingSubfacesQueue);
	vector<DartHandle> cavity[2];
	unsigned int missingSubfaceId;

	vector<Triangle> cdtFacetList;
	vector<DartHandle> lcc3CellsToBeRemoved;

	createEquivalentTetrahedralization(); 

	while (missingSubfacesQueue.size() != 0)
	{
		missingSubfaceId = missingSubfacesQueue.back();
		missingSubfacesQueue.pop_back();
		formCavity(cavity, missingSubfaceId, lcc3CellsToBeRemoved);
		
		for (unsigned int cavityId = 0; cavityId < 2; cavityId++)
		{
			cavityRetetrahedralization(cavity[cavityId], lcc3CellsToBeRemoved); 
			
		}
		
	}

	return;	
}		 
	


/////////////////////////////////////////////// Facet recovery ends /////////////////////////////////////////////////////////////
	

// main procedure
int main()
{
	readPLCInput();
	computeDelaunayTetrahedralization();
	recoverConstraintSegments();
	//removeLocalDegeneracies();
	//recoverConstraintFaces();

	return 0;
}

