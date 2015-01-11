#include <iostream>
#include <vector>
#include<math.h>

#include <CGAL/Object.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Circle_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Spherical_kernel_intersections.h>
#include "rply/rply.h"

#define Pi 22.0/7.0
#define INVALID_VALUE -1.0f // used in context of distances 
using namespace std;
using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;
typedef Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef Triangulation_data_structure_3<Vb> Tds;
typedef Delaunay_triangulation_3<K, Tds, Fast_location> Delaunay;
typedef Delaunay::Point Point;
typedef Sphere_3<K> Sphere;
typedef Circle_3<K> Circle;

typedef Exact_spherical_kernel_3 SK;
typedef Line_arc_3<SK> SphericalSegment;
typedef Sphere_3<SK> SphericalSphere;
typedef Segment_3<SK> CGALSegment; 
typedef Point_3<SK> SphericalPoint;

typedef Delaunay::Vertex_handle Vertex_handle;
/*
 * Input  : plcVertices, plcSegments, plcFaces
 * Output : cdtTeterahedralMesh (collection of tetrahedrons), cdtVertices(=plcVertices)
 */

class Segment
{
	public:
		unsigned int pointIds[2]; // index into the plcVertices vector
};

class Triangle
{
	public:
		unsigned int pointIds[3];
		
};

class Tetrahedron
{
	public:
		unsigned int pointIds[4];
		
};


vector<pair<Point, unsigned int> > plcVertices;
vector<Segment> plcSegments;
vector<Triangle> plcFaces;

vector <Tetrahedron> cdtTets; // contains ids of vertices making tetrahedron

Delaunay DT;

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
		plcVertices.push_back(make_pair(Point(tempPoint[0],tempPoint[1],tempPoint[2]), pointCount++));
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
			plcFaces.push_back(tempFace);
	                
			// extract corresponding segments as well
			Segment tempSegments[3]; // each face will give 3 segments
			
			tempSegments[0].pointIds[0] = tempFace.pointIds[0];
			tempSegments[0].pointIds[1] = tempFace.pointIds[1];
			
			tempSegments[1].pointIds[0] = tempFace.pointIds[1];
			tempSegments[1].pointIds[1] = tempFace.pointIds[2];
		
			tempSegments[2].pointIds[0] = tempFace.pointIds[2];
			tempSegments[2].pointIds[1] = tempFace.pointIds[3];


			plcSegments.push_back(tempSegments[0]);
			plcSegments.push_back(tempSegments[1]);		
			plcSegments.push_back(tempSegments[2]);
		
			break;
        	default: 
                	break;
        } 
    
	return 1;
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

	cout << "Number of vertices:" << plcVertices.size() << "\n";
	cout << "Number of faces:" << plcFaces.size() << "\n";
	cout << "Number of segments:" << plcSegments.size();
}

//////////////////////////////////////////////////// Done with reading input PLC ////////////////////////////////////////////////
// computes delaunay tetrahedralization
void computeDelaunayTetrahedralization()
{	
	// Add vertex id with each vertex which will be used for relating these vertices with plcSegments, plcFaces 
	DT.insert(plcVertices.begin(), plcVertices.end());

	cout << "\nDelaunay tetrahedralization computed!!\n";

	cout << "Number of vertices in Delaunay tetrahedralization:" << DT.number_of_vertices() << "\n";
	cout << "Number of tetrahedrons in Delaunay tetrahedralization:" << DT.number_of_cells() << "\n";

}



///////////////////////////////////////////// Segment recovery ///////////////////////////////////////////////////////////////////


void formMissingSegmentsQueue(vector<unsigned int> &missingSegmentQueue)
{
	// collect pointers to all segments in plcSegments which are not in DT	
	// while inserting points in DT set a id for each point

	missingSegmentQueue.clear();


	bool  segmentFound;

	for (unsigned int m = 0; m < plcSegments.size(); m++)
	{
		segmentFound = false;

		for (Delaunay::Finite_cells_iterator cit = DT.finite_cells_begin(); cit != DT.finite_cells_end(); cit++)
		{
			// since a tet is fully connected I only need to find if there is a tet containing both vertices			
			for (unsigned int i = 0; i < 4; i++)
			{
				if ((cit->vertex(i))->info() == plcSegments[m].pointIds[0])
					for (unsigned int j = 0; j != i && j < 4; j++)
					{	if ((cit->vertex(j))->info() == plcSegments[m].pointIds[1])
						{
							segmentFound = true;
							break;	// segment found!!
						}
					}
				if (segmentFound)
					break;
			}
		
			if (segmentFound)
				break;
		}
	
		if (!segmentFound)
			missingSegmentQueue.push_back(m);
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

void computeReferencePoint(Point *refPoint, unsigned int missingSegmentId)
{
	// compute a reference point p such that:
	// p encroaches missingSegment	
	// p has the largest circumradius of smallest circumsphere out of all encroaching vertices
	
	
	Point &A = plcVertices[plcSegments[missingSegmentId].pointIds[0]].first;
	Point &B = plcVertices[plcSegments[missingSegmentId].pointIds[1]].first;

	float missingSegmentLength = sqrt(pow((A.x() - B.x()), 2) + pow((A.y() - B.y()), 2) + pow((A.z() - B.z()), 2));
	
	float sphereRadius = (missingSegmentLength / 2.0);
	Point sphereCenter = Point((A.x() + B.x()) / 2.0,  (A.y() + B.y()) / 2.0, (A.z() + B.z()) / 2.0);

	Sphere smallestCircumsphere(sphereCenter, pow(sphereRadius, 2));


	float encroachingCandidateDistance;
	vector<float> circumradiusMap; // value is computed only for those points which are encroaching
	
	for (unsigned int n = 0; n < plcVertices.size(); n++)
	{
		if (n != plcSegments[missingSegmentId].pointIds[0] && n != plcSegments[missingSegmentId].pointIds[1])
		{
			
			encroachingCandidateDistance = sqrt(pow(sphereCenter.x() - (plcVertices[n].first).x(), 2) + pow(sphereCenter.y() - (plcVertices[n].first).y(), 2) + pow(sphereCenter.z() - (plcVertices[n].first).z(), 2));
			if (encroachingCandidateDistance <= sqrt(smallestCircumsphere.squared_radius()) / 2.0) // encroaching vertex
				circumradiusMap.push_back(computeCircumradius(A, B, plcVertices[n].first));
			else // not encroaching
				circumradiusMap.push_back(INVALID_VALUE);
		
		}
		else
			continue; // skip
	}

	// select the one with maximum circumradius.
	float maxCircumradius = INVALID_VALUE;
	
	for (unsigned int i = 0; i < plcVertices.size(); i++)
		if (circumradiusMap[i] != INVALID_VALUE)
			if (maxCircumradius < circumradiusMap[i])
			{
				maxCircumradius = circumradiusMap[i];
				refPoint = &(plcVertices[i].first);
			}	
	return;
}

float dotProduct (unsigned segment1Id, unsigned int segment2Id)
{
	Point segment1Vertex[2];
	Point segment2Vertex[2];

	for (unsigned int i = 0; i < 2; i++)
	{
		segment1Vertex[i] = plcVertices[plcSegments[segment1Id].pointIds[i]].first;
		segment2Vertex[i] = plcVertices[plcSegments[segment2Id].pointIds[i]].first;
	}

	Point vector1 = Point(segment1Vertex[0].x() - segment1Vertex[1].x(), segment1Vertex[0].y() - segment1Vertex[1].y(), segment1Vertex[0].z() - segment1Vertex[1].z());
	Point vector2 = Point(segment2Vertex[0].x() - segment2Vertex[1].x(), segment2Vertex[0].y() - segment2Vertex[1].y(), segment2Vertex[0].z() - segment2Vertex[1].z());


	float v1Dotv2 = vector1.x() * vector2.x() + vector1.y() * vector2.y() + vector1.z() * vector2.z();

	return v1Dotv2;

}

float vectorMagnitude(unsigned int inputSegmentId)
{
	Point segmentVertices[2];

	for (unsigned int i = 0; i < 2; i++)
		segmentVertices[i] = plcVertices[plcSegments[inputSegmentId].pointIds[i]].first;

	Point vector = Point(segmentVertices[0].x() - segmentVertices[1].x(), segmentVertices[0].y() - segmentVertices[1].y(), segmentVertices[0].z() - segmentVertices[1].z());

	float vectorMagnitude = sqrtf(powf(vector.x(), 2.0) + powf(vector.y(), 2.0) + powf(vector.z(), 2.0));

	return vectorMagnitude;
}


float computeAngleBetweenSegments(unsigned int segment1Id, unsigned int segment2Id)
{
	float angle = acosf(dotProduct(segment1Id, segment2Id) / (vectorMagnitude(segment1Id) * vectorMagnitude(segment2Id)));

	return (angle * 180.0f / Pi); // conversion to degrees
}

bool isVertexAcute(unsigned int inputPointId)
{
	// Determine segment-pair(involving A) 
	vector<unsigned int> incidentOnInputPoint;

	for (unsigned int n = 0; n < plcSegments.size(); n++) 
	{
		if (plcSegments[n].pointIds[0] == inputPointId || plcSegments[n].pointIds[1] == inputPointId)
			incidentOnInputPoint.push_back(n);
	}

	// Compute angle between all possible pairs(NAIVE SOLUTION)
	for (vector<unsigned int>::iterator segIter1 = incidentOnInputPoint.begin(); segIter1 != incidentOnInputPoint.end(); segIter1++) 
		for (vector<unsigned int>::iterator segIter2 = incidentOnInputPoint.begin(); *segIter1 != *segIter2 && segIter2 != incidentOnInputPoint.end(); segIter2++)
			if (computeAngleBetweenSegments(*segIter1, *segIter2) < 90.0f)
				return true;

	return false; 
}


unsigned int determineSegmentType(unsigned int missingSegmentId)
{

	// if both endpoints of the segment are acute, type 1
	// if only one endpoint acute, type 2	
	Point &A = plcVertices[plcSegments[missingSegmentId].pointIds[0]].first;
	Point &B = plcVertices[plcSegments[missingSegmentId].pointIds[1]].first;

	bool vertexAIsAcute = isVertexAcute(plcSegments[missingSegmentId].pointIds[0]);
	bool vertexBIsAcute = isVertexAcute(plcSegments[missingSegmentId].pointIds[1]);

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



bool containsSegment(unsigned int faceId, unsigned int segmentId)
{
	for (unsigned int i = 0; i < 3; i++)
		if (plcSegments[segmentId].pointIds[0] == plcFaces[faceId].pointIds[i])
			for (unsigned int n = 0; i !=n && n < 3; n++)
				if (plcSegments[segmentId].pointIds[2] == plcFaces[faceId].pointIds[n])
					return true;
	return false;
}


void updatePLCAndDT(Point v, unsigned int missingSegmentId)
{
	// vertex
	
	plcVertices.push_back(make_pair(v, plcVertices.size()));
	
	// segments
	Segment AN, NB;

	
	AN.pointIds[0] = plcSegments[missingSegmentId].pointIds[0];
	AN.pointIds[1] = plcVertices.size() - 1; // since new vertex has been inserted at the end
	NB.pointIds[0] = plcVertices.size() - 1;
	NB.pointIds[1] = plcSegments[missingSegmentId].pointIds[1];
	
	plcSegments.push_back(AN);
	plcSegments.push_back(NB);

	plcSegments.erase(plcSegments.begin() + missingSegmentId);

	// faces
	// Find out the faces sharing that segment
	 	//For each face sharing this segment edge,
			// Partition the face into 2 triangles
	unsigned int v1, v2, v3;
	Triangle newFace1, newFace2;

	for (unsigned int d = 0; d < plcFaces.size(); d++)
		if (containsSegment(d, missingSegmentId))
		{
			newFace1.pointIds[0] = plcFaces[d].pointIds[0];
			newFace1.pointIds[1] = plcFaces[d].pointIds[1];
			newFace1.pointIds[2] = plcVertices.size() - 1;
			
			newFace2.pointIds[0] = plcFaces[d].pointIds[0];
			newFace2.pointIds[1] = plcFaces[d].pointIds[2];
			newFace2.pointIds[2] = plcVertices.size() - 1;
	
			plcFaces.push_back(newFace1);
			plcFaces.push_back(newFace2);
			plcFaces.erase(plcFaces.begin() + d);		
		}


	// DT
	computeDelaunayTetrahedralization();


}


void splitMissingSegment(unsigned int missingSegmentId)
{

	// Apply rules to split the segments into strongly Delaunay subsegments
	Point vb, refPoint;
	Point sphereCenter;
	float sphereRadius;

	unsigned int segmentType;
	segmentType = determineSegmentType(missingSegmentId);

	computeReferencePoint(&refPoint, missingSegmentId);
	
	unsigned int A = plcSegments[missingSegmentId].pointIds[0];
	unsigned int B = plcSegments[missingSegmentId].pointIds[1];
	
	float AP, PB, AB;

	Point v;

	if (segmentType == 1)
	{
		Point &Areference = plcVertices[A].first;
		Point &Breference = plcVertices[B].first;
		AP = computeSegmentLength(Areference, refPoint);
		AB = computeSegmentLength(Areference, Breference);
		
		if (AP < 0.5f * AB)
		{
			sphereCenter = Areference;
			sphereRadius = AP;
		}

		else if (PB < 0.5 * AB)
		{
			sphereCenter = Breference;
			sphereRadius = PB;
		}

		else
		{
			sphereCenter = Areference;
			sphereRadius = 0.5 * AB; 
		}	

		SphericalPoint sphericalSphereCenter = SphericalPoint(sphereCenter.x(), sphereCenter.y(), sphereCenter.z());

		SphericalSphere s = SphericalSphere(sphericalSphereCenter, pow(sphereRadius, 2));

		SphericalPoint p1 = SphericalPoint(plcVertices[plcSegments[missingSegmentId].pointIds[0]].first.x(), plcVertices[plcSegments[missingSegmentId].pointIds[0]].first.y(), plcVertices[plcSegments[missingSegmentId].pointIds[0]].first.z());

		SphericalPoint p2 = SphericalPoint(plcVertices[plcSegments[missingSegmentId].pointIds[1]].first.x(), plcVertices[plcSegments[missingSegmentId].pointIds[1]].first.y(), plcVertices[plcSegments[missingSegmentId].pointIds[1]].first.z());

		CGALSegment seg(p1, p2);
		SphericalSegment lineArc(seg);
		vector<Object> intersections;
		CGAL::intersection(lineArc, s, back_inserter(intersections));
		if (intersections.size() > 0)
		{	v = object_cast<Point>(intersections.back());
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
		unsigned int acuteParentId; 
		Segment ApB; 
		float vbLength;
		A = plcSegments[missingSegmentId].pointIds[0];
		B = plcSegments[missingSegmentId].pointIds[1];

		if (isVertexAcute(A) && !isVertexAcute(B))
			acuteParentId = plcSegments[missingSegmentId].pointIds[0];
	        else if (isVertexAcute(B) && !isVertexAcute(A))
		{
			acuteParentId = plcSegments[missingSegmentId].pointIds[1];
			// swap A <-> B 
			unsigned int t;
			t = A;
			A = B;
			B = t;
		}
		
		
		// Segment calculations
		
		ApB.pointIds[0] = acuteParentId;
		ApB.pointIds[1] = B;
		
	
		SphericalPoint p1(plcVertices[acuteParentId].first.x(), plcVertices[acuteParentId].first.y(), plcVertices[acuteParentId].first.z());
		SphericalPoint p2(plcVertices[plcSegments[missingSegmentId].pointIds[1]].first.x(), plcVertices[plcSegments[missingSegmentId].pointIds[1]].first.y(), plcVertices[plcSegments[missingSegmentId].pointIds[1]].first.z());
		
		CGALSegment seg(p1, p2);
		SphericalSegment lineArc(seg);
		
		/// Sphere calculations
		SphericalPoint acuteParent = SphericalPoint(plcVertices[acuteParentId].first.x(), plcVertices[acuteParentId].first.y(), plcVertices[acuteParentId].first.z());
		
		Point acuteParentLinearField(plcVertices[acuteParentId].first.x(), plcVertices[acuteParentId].first.y(), plcVertices[acuteParentId].first.z());
	
		float ApRefPointLength = computeSegmentLength(acuteParentLinearField, refPoint);

		SphericalSphere s(acuteParent, pow(ApRefPointLength, 2));

		vector<Object> intersections;

		CGAL::intersection(lineArc, s, back_inserter(intersections));
	
		if (intersections.size() > 0)
		{	
			v = object_cast<Point>(intersections.back());
			intersections.pop_back();
		}

		else
		{
			cout << "Case 2: Sphere and segment do not intersect";
			cout << "\nSphere center:" << "(" << acuteParent.x() << "," << acuteParent.y() << "," << acuteParent.z() << ")" ;
			cout << "\nSphere radius:" << ApRefPointLength << "\n";
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
		static vector<unsigned int> acuteParentMap; // will store acute parent ID for each vertex
		unsigned int acuteParentId;
		// split the segment to create 2 type-2 segments:
		// Lets insert the new point at the mid of the segment
		Point newPoint;
		float x = (plcVertices[A].first.x() + plcVertices[B].first.x()) / 2.0;
		float y = (plcVertices[A].first.y() + plcVertices[B].first.y()) / 2.0;
		float z = (plcVertices[A].first.z() + plcVertices[B].first.z()) / 2.0;
		newPoint = Point(x, y, z); // must be collinear with A and B, we now have AX, XB
		
		v = newPoint;
	}

	// udpdate plc and DT
	
	return;
}


void recoverConstraintSegments()
{
	// I/P: plc1, DT1
	// O/P: plc2, DT2
	
	vector<unsigned int> missingSegmentQueue;// stores index of missing segments into plcSegments	
	unsigned int missingSegment;

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
	
	cout << "\nIn the loop" << i << "number of times" << "\n";
	
	return;
}

/////////////////////////////////////////////// Local Degeneracy Removal begin ///////////////////////////////////////////////////

class DegenerateVertexSetCandidate
{
	public:
		unsigned int pointIds[5];
};


// returns true if  given vertices are co-spherical
bool areCospherical(DegenerateVertexSetCandidate degenSet)
{	
	Point p[5];
	
	for (unsigned int i = 0; i < 5; i++)	
		p[i] = plcVertices[degenSet.pointIds[i]].first;

	if (CGAL::side_of_bounded_sphere(p[0],p[1],p[2],p[3],p[4]) == CGAL::ON_BOUNDARY)
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
				Vertex_handle vh = DT.mirror_vertex(cellIter, j);	
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
	
/*	for (unsigned int n = 0; n < 5; n++)
		if (isVertexPerturbable(degenCandidate.pointIds[n]))
			if (isVertexSegmentSafePerturbable())
			{
				removable = true;
				break;
			}
			else
				continue;

*/

	return removable;
}

void perturbRemove(unsigned int pointId, vector<DegenerateVertexSetCandidate>& localDegeneracySet)
{
	// DO NOTHING FOR NOW
}


bool areAffinelyIndependent(DegenerateVertexSetCandidate degenSet)
{

	// test each 4-tuple for coplanarity
	// If they are coplanar then return false
	// else return true
	Point v1 = plcVertices[degenSet.pointIds[0]].first;
	Point v2 = plcVertices[degenSet.pointIds[1]].first;
	Point v3 = plcVertices[degenSet.pointIds[2]].first;
	Point v4 = plcVertices[degenSet.pointIds[3]].first;
	Point v5 = plcVertices[degenSet.pointIds[4]].first;

	if (coplanar(v1, v2, v3, v4) || coplanar(v2, v3, v4, v5) || coplanar(v3, v4, v5, v1) || coplanar(v4, v5, v1, v2) || coplanar(v5, v1, v2, v3))
		return false;

	return true; // because no 4-tuple of points are coplanar
}

void computeBreakPoint(Point &vb, unsigned int pointId, vector<DegenerateVertexSetCandidate> &localDegeneracySet)
{
	// if vertices of localDeneracySet[n] are affinely independent
		// Let S be there common sphere, then vb is inside s
	// if vertices of localDegeneracySet[n] aren't affinely independent(4 of them are coplanar)
		// Let C be the common circle of those coplanar vertices, then vb lines inside C
	
	Point v1 = plcVertices[localDegeneracySet[pointId].pointIds[0]].first;
	Point v2 = plcVertices[localDegeneracySet[pointId].pointIds[1]].first;
	Point v3 = plcVertices[localDegeneracySet[pointId].pointIds[2]].first;
	Point v4 = plcVertices[localDegeneracySet[pointId].pointIds[3]].first;
	Point v5 = plcVertices[localDegeneracySet[pointId].pointIds[4]].first;

	if (areAffinelyIndependent(localDegeneracySet[pointId]))
	{
		// find common sphere
		// if I define a sphere using only 4 of these vertices then that will be shared by 5 as well:
	
		Sphere s(v1, v2, v3, v4);
		vb = s.center();
	}	

	else // 4 points are coplanar
	{
		// find common circle
		// find which 4 points are coplanar
		// find corresponding circle
		Circle c;

		if (coplanar(v1, v2, v3, v4))	
			c = Circle(v1, v2, v3); // v4 will also lie on this circle, since these vertices are co-spherical and coplanar.

		if (coplanar(v2, v3, v4, v5))	
			c = Circle(v2, v3, v4);

		if (coplanar(v3, v4, v5, v1))	
			c = Circle(v3, v4, v5);

		if (coplanar(v4, v5, v1, v2))	
			c = Circle(v4, v5, v1);

		if (coplanar(v5, v1, v2, v3))	
			c = Circle(v5, v1, v2);

		vb = c.center(); 
	}
}

bool isVertexEncroachingSegment(Point vb, unsigned int segmentId)
{
	// compute diametric sphere and check if:
		// radius < distance of center from vb
			// if yes, return true
			// else return false
	float x = (plcVertices[plcSegments[segmentId].pointIds[0]].first.x() + plcVertices[plcSegments[segmentId].pointIds[1]].first.x()) / 2.0;
	float y = (plcVertices[plcSegments[segmentId].pointIds[0]].first.y() + plcVertices[plcSegments[segmentId].pointIds[1]].first.y()) / 2.0;
	float z = (plcVertices[plcSegments[segmentId].pointIds[0]].first.z() + plcVertices[plcSegments[segmentId].pointIds[1]].first.z()) / 2.0;

	Point sphereCenter(x, y, z);

	float sphereRadius = sqrt(pow(sphereCenter.x() - plcVertices[plcSegments[segmentId].pointIds[0]].first.x(), 2) + pow(sphereCenter.y() - plcVertices[plcSegments[segmentId].pointIds[0]].first.y(), 2) + pow(sphereCenter.z() - plcVertices[plcSegments[segmentId].pointIds[0]].first.z(), 2));

	float centerToBreakingPointDistance = sqrt(pow(sphereCenter.x() - vb.x(), 2) + pow(sphereCenter.y() - vb.y(), 2) + pow(sphereCenter.z() - vb.z(), 2));


	if (sphereRadius <= centerToBreakingPointDistance)
		return true;
	else
		return false;
		
}

bool isVertexEncroachingFace(Point vb, unsigned int faceId)
{
	// make a diametric sphere of face(or triangle in our case)
	// check if radius < distance of breaking point from center
	
	// radius & center of circle of 3 points = radius & center of their diametric sphere 
	
	Point v1(plcVertices[plcFaces[faceId].pointIds[0]].first);
	Point v2(plcVertices[plcFaces[faceId].pointIds[1]].first);
	Point v3(plcVertices[plcFaces[faceId].pointIds[2]].first);

	Circle c(v1, v2, v3);

	float centerToBreakingPointDistance = sqrt(pow((c.center()).x() - vb.x(), 2) + pow((c.center()).y() - vb.y(), 2) + pow((c.center()).z() - vb.z(), 2));

	if (sqrtf(c.squared_radius()) < centerToBreakingPointDistance)
		return true;

	return false;
}


bool isEncroachingPLC(Point vb)
{
	// encroaches any segment
	for (unsigned int n = 0; n < plcSegments.size(); n++)
		if (isVertexEncroachingSegment(vb, n))
			return true;

	// encroaches any face
	for (unsigned int k = 0; k < plcFaces.size(); k++)
		if (isVertexEncroachingFace(vb, k))
			return true;

	return false;
}

void boundaryProtection(Point vb)
{
	// for each encroached segment:
		// add its perturbed circumcenter to plcVerices
		// update PLC & DT 
	// for each encroached face:
		// compute circumcenter x
		// if x encroaches any segment, don't insert, goto to segment spliting using step 1 above
		// else, add x to PLC & DT 	
		
	// call Delaunay segment recovery to recovery all missing segments(new segments might have been created due to new vertex insertions)
	
	for (unsigned int n = 0; n < plcSegments.size(); n++)
	{
		if (isVertexEncroachingSegment(vb, n))
		{
			// circumcenter of segment 'n':
			float x = (plcVertices[plcSegments[n].pointIds[0]].first.x() + plcVertices[plcSegments[n].pointIds[1]].first.x()) / 2.0;
			float y = (plcVertices[plcSegments[n].pointIds[0]].first.y() + plcVertices[plcSegments[n].pointIds[1]].first.y()) / 2.0;
			float z = (plcVertices[plcSegments[n].pointIds[0]].first.z() + plcVertices[plcSegments[n].pointIds[1]].first.z()) / 2.0;

			Point sphereCenter(x, y, z);

			updatePLCAndDT(sphereCenter, n);
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
				Point vb;	
				computeBreakPoint(vb, n, localDegeneracySet);
				if (isEncroachingPLC(vb))
					boundaryProtection(vb);
				else
					plcVertices.push_back(make_pair(vb, plcVertices.size()));		
			}						
		}
	}
	cout << "\nLocal degeneracy removal completed";
}

///////////////////////////////////////////////////Local Degeneracy Removal Ends//////////////////////////////////////////////////////



/////////////////////////////////////////////// Facet recovery starts ///////////////////////////////////////////////////////////

void formMissingSubfaceQueue(vector<unsigned int> &missingSubfacesQueue)
{
	// for all plcFaces check if those faces are already in DT 
	// If not, add them to missingSubfaceQueue
	// Else, continue
	Delaunay::Cell_handle ch;
	int i , j , k;
	
	Delaunay::Vertex v1, v2, v3;
	Vertex_handle vh1, vh2, vh3;

	for (unsigned int n = 0; n < plcFaces.size(); n++)
	{
		DT.is_vertex(plcVertices[plcFaces[n].pointIds[0]].first, vh1);
		DT.is_vertex(plcVertices[plcFaces[n].pointIds[1]].first, vh2);
		DT.is_vertex(plcVertices[plcFaces[n].pointIds[2]].first, vh3);

		if (DT.is_facet(vh1, vh2, vh3, ch, i, j, k))
			continue;
		else
			missingSubfacesQueue.push_back(n);
	}
}

void formCavity(vector<unsigned int>& cavity[2], unsigned int missingSubfaceId)
{
	// compute list of tets intersecting face number: missingSubfaceId
	// for each intersecting tet:
		// Remove the tet from DT 
		// Decrement counter of each facet(of neighboring tets) in contact of the removed tet by 1
		// For all removed tets, set counter of each of their facet = INVALID_VALUE
	// for all faces in DT:
		// Add all faces with counter value = 1 to cavity
		
	
	// Partition cavity into top & bottom cavities:
		// For each face in global cavity
			// Determine if it is on upper or lower side of missingSubface
				// Take any point on the face and test orientation of point wrt. missingSubface
					// Positive means associated facet belongs to upper cavity
					// Negative means associated facet belongs to lower cavity	 

	vector<Cell_iterator> intersectingTets;				
	for (Delaunay::Finite_cells_iterator cit = DT.finite_cells_begin(); cit != DT.finite_cells_end(); cit++)
	{
		if (areIntersecting(missingSubfaceId, cit))
			intersectingTets.push_back(cit);
		else
			continue;
	}		
					

}

// recovers the constraint faces
void recoverConstraintFaces()
{

	// input X2, D2
	// output CDT of X2
	
	// form a queue of missing subfaces
	// while (Q!=0)
	// 	remove an uncovered subface f from Q
	// 	for 2 cavities C1, C2 by formcavity procedure
	// 	for each cavity Ci:
	// 	call cavity retetrahedralization subroutine
	 			
        vector<unsigned int> missingSubfacesQueue;
	formMissingSubfaceQueue(missingSubfacesQueue);
	vector<Triangle> cavity[2];
	unsigned int missingSubfaceId;
	
	while (missingSubfacesQueue.size() != 0)
	{
		 missingSubfaceId = missingSubfacesQueue.back();
		 missingSubfacesQueue.pop_back();
		 formCavity(cavity, missingSubfaceId);
		
		 for (unsigned int cavityId = 0; cavityId < 2; cavityId++)
		 	cavityReterahedralization(cavity[cavityId], cdtTets); 
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
	removeLocalDegeneracies();
	recoverConstraintFaces();

	return 0;
}

