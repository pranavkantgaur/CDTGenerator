#include <string.h>
#include <iostream>
#include <unordered_set>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include "rply/rply.h"

#define Pi 22.0/7.0
#define INVALID_VALUE -1.0f // distances cannot be negative
using namespace std;
using namespace CGAL;


typedef Exact_predicates_inexact_constructions_kernel K;
typedef Point_3<K> Point;
typedef Delaunay_triangulation_3<K> Delaunay;
typedef Delaunay::Vertex_handle Vertex_handle;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PLC:
map <unsigned, Point> plcVertices; // mapping between vertex coordinates and corresponding unique id

class Segment
{
	unsigned int vertexIds[2]; // simply stores the endpoint ids
};

list<Segment> plcSegments;

class TriangleFace
{
	public:
		unsigned int vertexIds[3];
};

list <TriangleFace> plcFaces; // contains ids of vertices/points making the triangle

class TetrahedronCell
{
	public:
		unsigned int vertexIds[4];
};

map <unsigned, Point> cdtVertices;
list <TetrahedronCell> cdtTets; // contains ids of vertices making tetrahedron

Delaunay DT;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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
			plcFaces.push(tempFace);
	                
			// extract corresponding segments as well
			Segment tempSegments[3]; // each face will give 3 segments
			
			tempSegments[0].vertexIds[0] = tempFace.pointIds[0];
			tempSegments[0].vertexIds[1] = tempFace.pointIds[1];
			
			tempSegments[1].vertexIds[0] = tempFace.pointIds[1];
			tempSegments[1].vertexIds[1] = tempFace.pointIds[2];

			
			tempSegments[2].vertexIds[0] = tempFace.pointIds[2]
			tempSegments[2].vertexIds[1] = tempFace.pointIds[3];


			plcSegments.push(tempSegments[0]);
			plcSegments.push(tempSegments[1]);		
			plcSegments.push(tempSegments[2]);
		
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
//////////////////////////////////////////////////// Done with reading input PLC ////////////////////////////////////////////////
// computes delaunay tetrahedralization
void computeDelaunayTetrahedralization()
{	
	
	list <Point> tempPointList;
	unsigned int i = 0;

	for (map<unsigned, Point>::iterator pit = plcVertices.begin(); pit != plcVertices.end(); pit++, i++)
		tempPointList.push(plcVertices.find(i)->second);

    	DT.insert(tempPointList.begin(), tempPointList.end());

	cout << "\nInitial Delaunay tetrahedralization computed!!";
}
///////////////////////////////////////////// Segment recovery(ensures that all constraining segments are strongly Delaunay)/////
void formMissingSegmentsQueue(queue<Segment>missingSegmentQueue)
{
	// plcSegments contains the constraint segments
	// simply insert all segments from plcSegments to missingSegmentQueue
	missingSegmentQueue = plcSegments;

	return;
}

unsigned int computeCircumradius(Vertex A, Vertex B, Vertex encroachingCandidate)
{
	// computing circumradius for a triangle:
	
	float a, b, c;

	a = sqrt(pow((A.x - encroachingCandidate.x), 2)+ pow((A.y - encroachingCandidate.y), 2) + pow((A.z - encroachingCandidate.z), 2));
	b = sqrt(pow((B.x - encroachingCandidate.x), 2)+ pow((B.y - encroachingCandidate.y), 2) + pow((B.z - encroachingCandidate.z), 2));
	c = sqrt(pow((A.x - B.x), 2)+ pow((A.y - B.y), 2) + pow((A.z - B.z), 2));

	return circumradius = (a * b * c) / sqrt((a + b + c) * (b + c - a) * (c + a - b) * (a + b - c));
}

void computeReferencePoint(Vertex *refPoint, Segment *missingSegment)
{
	// compute a reference point p such that:
	// p encroaches missingSegment	
	// p has the largest circumradius of smallest circumsphere out of all encroaching vertices
	
	Vertex A, B;
	getVertexbyId(missingSegment->vertexId[0], A);
	getVertexbyId(missingSegment->vertexId[1], B);

	float missingSegmentLength = sqrt(pow((A.x - B.x), 2) + pow((A.y - B.y), 2) + pow((A.z - B.z), 2));
	
	Sphere smallestCircumsphere;

	smallestCircumsphere.center.x = (A.x + B.x) / 2;
	smallestCircumsphere.center.y = (A.y + B.y) / 2;
	smallestCircumsphere.center.z = (A.z + B.z) / 2;
	smallestCircumsphere.radius = (missingSegmentLength / 2.0);

	float encroachingCandidateDistance;
	float circumradiusMap[plcVertices.size() - 2]; // stores the circumradius value for all vertices for the triangle formed by missingsegment and vertex.
	for (unsigned int n = 0; n < plcVertices.size(); n++)
	{
		if (n != missingSegment->verexIds[0] && n != missingSegment->vertexIds[1])
			if (encroachingCandidateDistance <= smallestCircumsphere.radius) // encroaching vertex
				circumradiusMap[n] = computeCircumradius(A, B, plcVertices[n]);
			else // not encroaching
				circumradiusMap[i] = INVALID_VALUE;
		else
			continue; // skip
	}

	// select the one with maximum circumradius.
	float maxCircumradius = INVALID_VALUE;


	for (unsigned int i = 0; i < plcVertices.size(); i++)
		if (maxCircumradius < circumradius[i])
		{
			maxCircumradius = circumradius[i];
			refPoint = &plcVertices[i]; // at the end of iteration we will have reference point pointed to be refPoint
		}	
	return;
}

// returns vertex handle corresponding vertexId
Vertex& getVertexbyId(unsigned int vertexId)
{
	return plcVertices[vertexId];
}



float dotProduct (Segment segment1, Segment segment2)
{
	Vertex segment1Vertex[2];
	Vertex segment2Vertex[2];

	segment1Vertex[0] = getVertexbyId(segment1.vertexIds[0]);
	segment1Vertex[1] = getVertexbyId(segment1.vertexIds[1]);

	segment2Vertex[0] = getVertexbyId(segment2.vertexIds[0]);
	segment2Vertex[1] = getVertexbyId(segment2.vertexIds[1]);

	Vector vector1 = (segment1Vertex[0].x - segment1Vertex[1].x, segment1Vertex[0].y - segment1Vertex[1].y, segment1Vertex[0].z - segment1Vertex[1].z);
	Vector vector2 = (segment2Vertex[0].x - segment2Vertex[1].x, segment2Vertex[0].y - segment2Vertex[1].y, segment2Vertex[0].z - segment2Vertex[1].z);


	float v1Dotv2 = vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;

	return v1Dotv2;

	
}

float vectorMagnitude(Segment inputSegment)
{
	Vertex segmentVertices[2];
	segmentVertices[0] = getVertexbyId(inputSegment.vertexIds[0]);
	segmentVertices[1] = getVertexbyId(inputSegment.vertexIds[1]);

	Vector vector = (segmentVertices[0].x - segmentVertices[1].x, segmentVertices[0].y - segmentVertices[1].y, segmentVertices[0].z - segmentVertices[1].z);

	float sqrt(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2));
}


float computeAngleBetweenSegments(Segment segment1, Segment segment2)
{
	float angle = acosf(dotProduct(segment1, segment2) / (vectorMagnitude(segment1) * vectorMagnitude(segment2)));

	return (angle * 180.0f / Pi); // convertion to degrees
}

bool isVertexAcute(Vertex A)
{
	// Determine segment-pair(involving A) 
	list<Segment> incidentOnA;

	for (unsigned int n = 0; n < plcSegments.size(); n++) 
	{
		if (plcSegments[n].vertexIds[0] == A.id || plcSegments[n].vertexIds[1] == A.id)
			incidentOnA.push(plcSegments[n]);
	}

	// Compute angle between all possible pairs(NAIVE SOLUTION)
	for (list<Segment>::iterator segIter1 = incidentOnA.begin(); segIter1 != incidentOnA.end(); segIter1++) 
		for (list<Segment>::iterator segIter2 = incidentOnA.begin(); segIter2 != incidentOnA.end(); segIter2++)
			if ((segIter1 != segIter2) && (computeAngleBetweenSegments(segIter1, segIter2) < 90.0f))
				return true;

	return false; // statement is outside 'for' structure
}


unsigned int determineSegmentType(Segment *missingSegment)
{

	// if both endpoints of the segment are acute, type 1
	// if only one endpoint acute, type 2	
	Vertex A = getVertexbyId(missingSegment->vertexIds[0]);
	Vertex B = getVertexbyId(missingSegment->vertexIds[1]);

	bool vertexAIsAcute = isVertexAcute(A);
	bool vertexBIsAcute = isVertexAcute(B);

	if (vertexAIsAcute && vertexBIsAcute)	
		return 1;
	
	if (!vertexAIsAcute && !vertexBIsAcute)
		return 0; // invalid type

	else 
		return 2;
}

unsigned int findAcuteParent(unsigned int vertexId)
{
	if (isVertexAcute(vertexId))
		return vertexId;
	else
			

}

void splitMissingSegment(Segment *missingSegment)
{

	// Apply rules to split the segments into strongly Delaunay subsegments
	Vertex vb, refPoint;
	Sphere s;

	unsigned int segmentType;
	
	computeReferencePoint(refPoint);
	segmentType = determineSegmentType(missingSegment);

	Vertex A = getVertexbyId(missingSegment->vertexIds[0]);
	Vertex B = getVertexbyId(missingSegment->vertexIds[1]);
	
	float AP, PB, AB;

	Vertex v;

	if (segmentType == 1)
	{
		AP = computerSegmentLength(A, refPoint);
		AB = computerSegmentLength(A, B);
		
		if (AP < 0.5f * AB)
		{
			s.center = A;
			s.radius = AP;
		}

		else if (PB < 0.5 * AB)
		{
			s.center = B;
			s.radius = PB;
		}

		else
		{
			s.center = A;
			s.radius = 0.5 * AB; 
		}	

		v = CGAL::intersection(missingSegment, s);
	}

	else if (segmentType == 2)
	{
		// locate which vertex is not acute
		// lets determine which vertex is acute, if any
		// let vertex B be acute, then acute(B) = B, else acute(B) = parentVertex 
		unsigned int acuteParentId; 
		Segment ApB; 
		unsigned int vbLength;

	
		if (isVertexAcute(missingSegment->vertexIds[0]))
			acuteParentId = missingSegment->vertexIds[0];
	        else if (isVertexAcute(missingSegment->vertexIds[1]))
	       	       acuteParentId = missingSegment->vertexIds[1];
       		else
	        {
	        	if (A.acuteParentId != NULL)
		       	{
				acuteParentId = A.acuteParentId;
				ApB.vertexIds[1] = missingSegment->vertexIds[1]; // ie. vertex B taken as the other endpoint
		 		vbLength = computerSegmentLength(v, B);      	       	
		       	} 	       		

		       	else if (B.acuteParentId != NULL)
			{
				acuteParentId = B.acuteParentId;		
				ApB.vertexIds[1] = missingSegment->vertexIds[0]; // vertex A taken as the other endpoint
       	        		vbLength = computerSegmentLength(v, A);
			}
		}

		ApB.vertexIds[0] = acuteParentId;
	
		Vertex *acuteParent = getVertexbyId(acuteParentId);

		unsigned int ApRefPointLength = computerSegmentLength(acuteParent, refPoint);

		s.center = acuteParent;
		s.radius = ApRefPointLength;

		v = CGAL::intersection(ApB, s); 
		
		unsigned int vrefpointLength = computerSegmentLength(v, refPoint);

		if (vbLength < vrefpointLength) // v was rejected
		{
			s.center = acuteParent;
			unsigned int avLength = computerSegmentLength(A, v);
			if (vrefpointLength < 0.5 * avLength)
				{
			
					unsigned int acuteparentALength = computerSegmentLength(acuteParent, A);
					s.radius = acuteparentALength + avLength - vrefpointLength;	
				}
			else
				s.radius = acuteparentALength + 0.5 * avLength;
		
			
			v = CGAL::intersection(ApB, s);
		}
	}	
	
	// update PLC and Delaunay tetrahedralization
	plcVertices.push(v);
	computeDelaunayTetrahedralization();

	return;
}


void recoverConstraintSegments()
{
	// I/P: plcVertex1, plcFaces1, DT1
	// O/P: plcVertex2, plcFaces2, DT2
	
	formMissingSegmentsQueue(missingSegmentQueue);
 	Segment *missingSegment;

	while (missingSegmentQueue.size() != NULL)
	{
		missingSegment = missingSegmentQueue.pop();
		splitMissingSegment(missingSegment);
		updatePLCDT();
	}

	return;
}




/////////////////////////////////////////////// Local Degeneracy Removal begin ///////////////////////////////////////////////////


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
	DegenerateVertexSetCandidate degenerateSetCandidate;
	
	for (delaunay::Finite_cell_iterator cellIter = DT.finite_cells_begin(); cellIter != DT.finite_cells_end(); cellIter++)
	{
		for (unsigned int n = 0; n < 4; n++)
			degenerateSetCandidate.degenSetVertex[n] = cellIter->vertex(n);

		for (unsigned int j = 0; j < 4; j++)
			for (unsigned int k = 0; k < 4; k++)
				{
					if ((cellId->neighbor(j))->neighbor(k) == cellId)
						degenerateSetCandidate.degenSetVertex[4] = (cellIter->neighbor)->vertex(k);		
						
					if (areCospherical(degenerateSetCandidate))
					{
						if (localDegeneracySet.find(degenerateSetCandidate) != localDegeneracySet.end())
							localDegeneracySet.insert(degenerateSetCandidate);	
					}
				}
	}
}





// perturbation should be such that it does not make PLC inconsistent
bool isVertexPerturbable(Vertex)
{
	bool pertubable = false;
	
	// a vertex is perturbable iff its perturbation does not make PLC inconsitent
	// use simulation of simplicity approach to artifically perturb the vertex(if it is perturbable) to remove the degeneracy
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

///////////////////////////////////////////////////Local Degeneracy Removal Ends//////////////////////////////////////////////////////



/////////////////////////////////////////////////// Facet recovery starts ///////////////////////////////////////////////////////////

// IMPLEMENT UNORDERED_SET OF missingSubfacesQueue ????

void formMissingSubfaceQueue(unordered_set<> missingSubfacesQueue)
{
	
}

void formCavity(lcc c1, lcc c2, lcc aMissingSubface)
{

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
	 			
        unordered_set<face, key> missingSubfacesQueue;
	formMissingSubfaceQueue(missingSubfacesQueue);
	lcc aMissingSubface;
	lcc C[2]; // seems like we need to implement it using LCC


	while (missingSubfacesQueue.size() != 0)
	{
		 aMissingSubface = missingSubfacesQueue.pop();

		 formCavity(C[0], C[1], aMissingSubface);
		
		 for (unsigned int cavityId = 0; cavityId < 2; cavityId++)
		 	cavityReterahedralization(C[cavityId], cdtMesh); // 'cdtMesh' points to the final output mesh 
			 
	}

	return;	
}		 
	


//////////////////////////////////////////////////// Facet recovery ends /////////////////////////////////////////////////////////////
	

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

