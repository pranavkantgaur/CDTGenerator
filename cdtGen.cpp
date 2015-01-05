#include <iostream>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Sphere_3.h>
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


void formMissingSegmentsQueue(vector<unsigned int> missingSegmentQueue)
{
	// collect pointers to all segments in plcSegments which are not in DT	
	// while inserting points in DT set a id for each point
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

	a = sqrt(pow((A.x() - encroachingCandidate.x()), 2)+ pow((A.y() - encroachingCandidate.y()), 2) + pow((A.z() - encroachingCandidate.z()), 2));
	b = sqrt(pow((B.x() - encroachingCandidate.x()), 2)+ pow((B.y() - encroachingCandidate.y()), 2) + pow((B.z() - encroachingCandidate.z()), 2));
	c = sqrt(pow((A.x() - B.x()), 2)+ pow((A.y() - B.y()), 2) + pow((A.z() - B.z()), 2));

	return circumradius = (a * b * c) / sqrt((a + b + c) * (b + c - a) * (c + a - b) * (a + b - c));
}

void computeReferencePoint(Point *refPoint, Segment *missingSegment)
{
	// compute a reference point p such that:
	// p encroaches missingSegment	
	// p has the largest circumradius of smallest circumsphere out of all encroaching vertices
	
	Point &A, &B;
	A = plcVertices[missingSegment->pointIds[0]];
	B = plcVertices[missingSegment->pointIds[1]];

	float missingSegmentLength = sqrt(pow((A.x() - B.x()), 2) + pow((A.y() - B.y()), 2) + pow((A.z() - B.z()), 2));
	
	float sphereRadius = (missingSegmentLength / 2.0);
	Point sphereCenter = ((A.x() + B.x()) / 2,  (A.y() + B.y()) / 2, (A.z() + B.z()) / 2);

	Sphere smallestCircumsphere(sphereCenter, pow(sphereRadius, 2));


	float encroachingCandidateDistance;
	vector<float> circumradiusMap; // value is computed only for those points which are encroaching
	
	for (unsigned int n = 0; n < plcVertices.size(); n++)
	{
		if (n != missingSegment->pointIds[0] && n != missingSegment->pointIds[1])
		{
			
			encroachingCandidateDistance = sqrt(pow(sphereCenter.x() - plcVertices[n].x(), 2) + pow(sphereCenter.y() - plcVertices[n].y(), 2) + pow(sphereCenter.z() - plcVertices[n].z(), 2));
			if (encroachingCandidateDistance <= sqrt(smallestCircumsphere.squared_radius()) / 2.0) // encroaching vertex
				circumradiusMap[n] = computeCircumradius(A, B, plcVertices[n]);
			else // not encroaching
				circumradiusMap[i] = INVALID_VALUE;
		
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
				refPoint = &plcVertices[i];
			}	
	return;
}

float dotProduct (Segment &segment1, Segment &segment2)
{
	Point &segment1Vertex[2];
	Point &segment2Vertex[2];

	for (unsigned int i = 0; i < 2; i++)
	{
		segment1Vertex[i] = plcVertices[segment1.pointIds[i]];
		segment2Vertex[i] = plcVertices[segment2.pointIds[i]];
	}

	Point  vector1 = (segment1Vertex[0].x() - segment1Vertex[1].x(), segment1Vertex[0].y() - segment1Vertex[1].y(), segment1Vertex[0].z() - segment1Vertex[1].z());
	Point vector2 = (segment2Vertex[0].x() - segment2Vertex[1].x(), segment2Vertex[0].y() - segment2Vertex[1].y(), segment2Vertex[0].z() - segment2Vertex[1].z());


	float v1Dotv2 = vector1.x() * vector2.x() + vector1.y() * vector2.y() + vector1.z() * vector2.z();

	return v1Dotv2;

}

float vectorMagnitude(Segment &inputSegment)
{
	Vertex &segmentVertices[2];

	for (unsigned int i = 0; i < 2; i++)
		segmentVertices[i] = plcVertices[inputSegment.pointIds[i]];

	Point vector = (segmentVertices[0].x() - segmentVertices[1].x(), segmentVertices[0].y() - segmentVertices[1].y(), segmentVertices[0].z() - segmentVertices[1].z());

	float vectorMagnitude = sqrt(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2));

	return vectorMagnitude;
}


float computeAngleBetweenSegments(Segment *segment1, Segment *segment2)
{
	float angle = acosf(dotProduct(*segment1, *segment2) / (vectorMagnitude(*segment1) * vectorMagnitude(*segment2)));

	return (angle * 180.0f / Pi); // conversion to degrees
}

bool isVertexAcute(unsigned int inputPointId)
{
	// Determine segment-pair(involving A) 
	vector<Segment*> incidentOnInputPoint;

	for (unsigned int n = 0; n < plcSegments.size(); n++) 
	{
		if (plcSegments[n].pointIds[0] == inputPointId || plcSegments[n].pointIds[1] == inputPointId)
			incidentOnInputPoint.push_back(&plcSegments[n]);
	}

	// Compute angle between all possible pairs(NAIVE SOLUTION)
	for (vector<Segment*>::iterator segIter1 = incidentOnInputPoint.begin(); segIter1 != incidentOnInputPoint.end(); segIter1++) 
		for (vector<Segment*>::iterator segIter2 = incidentOnInputPoint.begin(); segIter2 != incidentOnInputPoint.end(); segIter2++)
			if ((segIter1 != segIter2) && (computeAngleBetweenSegments(segIter1, segIter2) < 90.0f))
				return true;

	return false; 
}


unsigned int determineSegmentType(Segment *missingSegment)
{

	// if both endpoints of the segment are acute, type 1
	// if only one endpoint acute, type 2	
	Point &A = plcVertices[missingSegment->pointIds[0]];
	Point &B = plcVertices[missingSegment->pointIds[1]];

	bool vertexAIsAcute = isVertexAcute(missingSegment->pointIds[0]);
	bool vertexBIsAcute = isVertexAcute(missingSegment->pointIds[1]);

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


void splitMissingSegment(unsigned int missingSegmentId)
{

	// Apply rules to split the segments into strongly Delaunay subsegments
	Point vb, refPoint;
	Point sphereCenter;
	float sphereRadius;

	unsigned int segmentType;
	segmentType = determineSegmentType(missingSegment);

	computeReferencePoint(&refPoint, missingSegment);
	
	Point &A = plcVertices[plcSegments[missingSegmentId].pointIds[0]];
	Point &B = plcVertices[plcSegments[missingSegmentId].pointIds[1]];
	
	float AP, PB, AB;

	Point v;

	if (segmentType == 1)
	{
		AP = computeSegmentLength(A, refPoint);
		AB = computeSegmentLength(A, B);
		
		if (AP < 0.5f * AB)
		{
			sphereCenter(A);
			sphereRadius = AP;
		}

		else if (PB < 0.5 * AB)
		{
			sphereCenter(B);
			sphereRadius = PB;
		}

		else
		{
			sphereCenter(A);
			sphereRadius = 0.5 * AB; 
		}	

		Sphere s(sphereCenter, pow(sphereRadius, 2));

		v = CGAL::intersection(plcSegments[missingSegmentId], s);
	}

	else if (segmentType == 2)
	{
		unsigned int acuteParentId; 
		Segment ApB; 
		float vbLength;


		if (isVertexAcute(A) && !isVertexAcute(B))
			acuteParentId = plcSegments[missingSegmentId].pointIds[0];
	        else if (isVertexAcute(B) && !isVertexAcute(A))
		{
			acuteParentId = plcSegments[missingSegmentId].pointIds[1];
			// swap A <-> B 
			Point t;
			t = A;
			A = B;
			B = t;
		}
		
		ApB.pointIds[0] = acuteParentId;
		ApB.pointIds[1] = B;
	
		Point &acuteParent = plcVertices[acuteParentId];

		float ApRefPointLength = computerSegmentLength(acuteParent, refPoint);

		Sphere s(acuteParent, pow(ApRefPointLength, 2));

		v = CGAL::intersection(ApB, s); 
		
		float vrefpointLength = computerSegmentLength(v, refPoint);

		if (vbLength < vrefpointLength) // v was rejected
		{
			sphereCenter = acuteParent;
			unsigned int avLength = computerSegmentLength(A, v);
			if (vrefpointLength < 0.5 * avLength)
				{
			
					float acuteparentALength = computerSegmentLength(acuteParent, A);
					sphereRadius = acuteparentALength + avLength - vrefpointLength;	
				}
			else
				sphereRadius = acuteparentALength + 0.5 * avLength;
		        
			s = Sphere(sphereCenter, pow(sphereRadius, 2));
			v = CGAL::intersection(ApB, s);
		}
	}

	else if (segmentType == 3) // meaning both A & B are acute
	{
		static vector<unsigned int> acuteParentMap; // will store acute parent ID for each vertex
		unsigned int acuteParentId;
		// split the segment to create 2 type-2 segments:
		// Lets insert the new point at the mid of the segment
		Point newPoint;
		x = (A.x() + B.x()) / 2.0;
		y = (A.y() + B.y()) / 2.0;
		z = (A.z() + B.z()) / 2.0;
		newPoint = Point(x, y, z); // must be collinear with A and B, we now have AX, XB
		
		v = newPoint;
	}

	// udpdate plc and DT
	
	// vertex
	plcVertices.push_back(v);
	
	// segments
	Segment AN, NB;

	
	AN.pointIds[0] = plcSegments[missingSegmentId].pointIds[0];
	AN.pointIds[1] = plcVertices.size() - 1; // since new vertex has been inserted at the end
	NB.pointIds[0] = plcVertices.size() - 1;
	NB.pointIds[1] = plcSegments[missingSegmentId].pointIds[1];
	
	plcSegments.push_back(AN);
	plcSegments.push_back(NB);

	for (unsigned int n = 0; n < plcSegments.size(); n++)
		if (plcSegments[n] == missingSegment)
		{
			plcSegments.erase(n); 
			break;
		}

	// faces
	// Find out the faces sharing that segment
	 	//For each face sharing this segment edge,
			// Partition the face into 2 triangles
			

	unsigned int oldFaceId;
	plcFaces.erase(oldFaceId);
	plcFaces.push_back(newFace1);
	plcFaces.push_back(newFace2);

	// DT
	computeDelaunayTetrahedralization();

	return;
}


void recoverConstraintSegments()
{
	// I/P: plc1, DT1
	// O/P: plc2, DT2
	
	vector<unsigned int> missingSegmentQueue;// stores index of missing segments into plcSegments	
	unsigned int missingSegment;

	formMissingSegmentsQueue(missingSegmentQueue);

	while (missingSegmentQueue.size() != 0)
	{
		missingSegment = missingSegmentQueue.pop_back();
		splitMissingSegment(missingSegment);
		formMissingSegmentsQueue(missingSegmentQueue);
	}

	return;
}


/*

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
	
*/

//////////////////////////////////////////////////// Facet recovery ends /////////////////////////////////////////////////////////////
	

// main procedure
int main()
{
	readPLCInput();
	computeDelaunayTetrahedralization();
	recoverConstraintSegments();
/*	
	removeLocalDegeneracies();
	recoverConstraintFaces();
*/
	return 0;
}

