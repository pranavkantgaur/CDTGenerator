#include "cdtGen.h"

LCC plc; /*!< Input piecewise linear cell complex representing the input*/
Delaunay DT; /*< Intermidiate structure used for storing Delaunay tetrahedralization*/
vector <CGALPoint> plcVertexVector; /*> Used for initializing plc*/

/*! \class Triangle
    \brief Represents triangle 
  
    Triangle class stores indices of vertices. 
*/
class Triangle
{
	public:
		size_t pointIds[3]; /*< Indices of points */
};


vector <Triangle> plcFaceVector; /*> Used for initializing plc*/

LCC cdtMesh; /*!< Output mesh */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


static float tempPoint[3]; /*< scratchpad for storing point */
static int pointCount = 0; /*< count for number of points */
unsigned int dimensionId = 0; /*< index into x, y, z component of a point */
static unsigned pointId = 0; /*< unique index of a point */

/*! \fn static int vertex_cb(p_ply_argument argument)
    \brief Callback for reading vertex from PLY file	
    
    \param argument represents PLY file
 */
static int vertex_cb(p_ply_argument argument) 
{	
	long eol;
	ply_get_argument_user_data(argument, NULL, &eol);
	tempPoint[dimensionId++] = ply_get_argument_value(argument);
	
	if (eol)
	{
		plcVertexVector.push_back(CGALPoint(tempPoint[0],tempPoint[1],tempPoint[2]));
		dimensionId = 0;
	}
	
	return 1;
}

/*! \fn static int face_cb(p_ply_argument argument)
    \brief Reads face information from PLY file

   \param argument Represents PLY file 
*/

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

/*! \fn bool areGeometricallySameSegments(DartHandle d1, DartHandle d2)
    \brief Tests whether segments represented by d1 and d2 are _geometrically_ same.
    \param [in] d1 DartHandle for first segment
    \param [in] d2 DartHandle for second segment
*/
bool areGeometricallySameSegments(DartHandle d1, DartHandle d2)
{
	if (plc.point(d1) == plc.point(plc.beta(d2, 1)))
		if (plc.point(plc.beta(d1, 1)) == plc.point(d2))
			return true;
	return false;
} 


/*! \fn void readPLCInput()
    \brief Reads PLC from input(only PLY supported currently) file

    Reads vertex and face information and initializes corresponding linear cell complex _plc_
 */
void readPLCInput()
{
	// read PLY file(assumed to contain the PLC)
	string fileName;

	cout << "\nPlease enter input filename:\t";
	cin >> fileName;
	
	p_ply inputPLY = ply_open(fileName.c_str(), NULL, 0, NULL);
    	
	if (!inputPLY) exit(0);

        if (!ply_read_header(inputPLY)) exit(0);

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
 
	int sewedMark = plc.get_new_mark();

	if (sewedMark == -1)
	{
		cout << "\nNo free mark available!!";
		exit(0);
	}
	
	CGALPoint trianglePoints[3];

	for (unsigned int n = 0, m = plcFaceVector.size(); n < m; n++)
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
					if (areGeometricallySameSegments(segIter1, segIter2)) // checks the geometry of segments
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


/*! \fn void writePLYOutput(LCC &lcc, string fileName)
    \brief Writes LCC to a PLY file

     Initializes PLY file header, vertices, faces from LCC and writes it to the location pointed to by _filename_

    \param [in] lcc Input linear cell complex
    \param [in] fileName Name of the putput file
 */
/*
 void writePLYOutput(LCC &lcc, string fileName)
{
	p_ply lccOutputPLY;

	if ((lccOutputPLY = ply_create(fileName.c_str(), PLY_ASCII, NULL, 0, NULL)) == NULL)
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
*/
/*! \fn void computeDelaunayTetrahedralization()
    \brief Computes Delaunay tetrahedralization
 */
void computeDelaunayTetrahedralization()
{	
	vector <pair <CGALPoint, DartHandle> > lccVertexVector;

	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
		lccVertexVector.push_back(pair<CGALPoint, DartHandle>(plc.point(pIter), pIter));

	DT.insert(lccVertexVector.begin(), lccVertexVector.end());

	cout << "\nDelaunay tetrahedralization computed!!";
	cout << "\nNumber of vertices in Delaunay tetrahedralization:" << DT.number_of_vertices();
	cout << "\nNumber of tetrahedrons in Delaunay tetrahedralization:" << DT.number_of_cells();
}



/*! \fn void formMissingSegmentsQueue(vector<DartHandle> &missingSegmentQueue)
    \brief Collects missing constraint segments in a queue

    Determines which constraint segments are missing from Delaunay tetrahedralization of vertices of PLC

    \param missingSegmentQueue [Out] Contains the _DartHandle_ to the missing constraint segments after execution.
*/
void formMissingSegmentsQueue(vector<DartHandle> &missingSegmentQueue)
{
	missingSegmentQueue.clear();

	for (LCC::One_dart_per_cell_range<1>::iterator segmentIter = plc.one_dart_per_cell<1>().begin(), segmentIterEnd = plc.one_dart_per_cell<1>().end(); segmentIter != segmentIterEnd; segmentIter++)
	{
		CGALPoint p1 = plc.point(segmentIter);
		CGALPoint p2 = plc.point(plc.beta(segmentIter, 1));

		CellHandle c;
		int i, j;
		VertexHandle vh1, vh2;

		if (DT.is_vertex(p1, vh1))
			if (DT.is_vertex(p2, vh2))
				if (!DT.is_edge(vh1, vh2, c, i, j))
					missingSegmentQueue.push_back(segmentIter);
	}

	cout << "\nTotal number of missing constraint segments:" << missingSegmentQueue.size() << "\n";


	return;
}

/*! \fn unsigned int computeCircumradius(CGALPoint &A, CGALPoint &B, CGALPoint &encroachingCandidate)
    \brief Computes circumradius of circle defined by input points

    \param [in] A Endpoint 1 of constraint segment
    \param [in] B Endpoint 2 of constraint segment
    \param [in] encroachingCandidate Candidate for reference point
*/
unsigned int computeCircumradius(CGALPoint &A, CGALPoint &B, CGALPoint &encroachingCandidate)
{
	// computing circumradius of a triangle:
	
	float a, b, c;
	float circumradius;
	a = sqrt(pow((A.x() - encroachingCandidate.x()), 2)+ pow((A.y() - encroachingCandidate.y()), 2) + pow((A.z() - encroachingCandidate.z()), 2));
	b = sqrt(pow((B.x() - encroachingCandidate.x()), 2)+ pow((B.y() - encroachingCandidate.y()), 2) + pow((B.z() - encroachingCandidate.z()), 2));
	c = sqrt(pow((A.x() - B.x()), 2)+ pow((A.y() - B.y()), 2) + pow((A.z() - B.z()), 2));

	return circumradius = (a * b * c) / sqrt((a + b + c) * (b + c - a) * (c + a - b) * (a + b - c));
}


/*! \fn void computeReferencePoint(CGALPoint *refPoint, DartHandle missingSegmentHandle)
    \brief Computes the reference point for segment splitting

    Computes reference point which will be used for determing the location of _steiner point_ for the input missing constraint segment.
    \param [out] refPoint Pointer to the computed reference point.
    \param [in] missingSegmentHandle Dart handle of missing constraint segment.
*/
void computeReferencePoint(CGALPoint *refPoint, DartHandle missingSegmentHandle)
{

	CGALPoint &A = plc.point(missingSegmentHandle);
	CGALPoint &B = plc.point(plc.beta(missingSegmentHandle, 1));

	float missingSegmentLength = sqrt(pow((A.x() - B.x()), 2) + pow((A.y() - B.y()), 2) + pow((A.z() - B.z()), 2));
	float sphereRadius = (missingSegmentLength / 2.0);
	CGALPoint sphereCenter = CGALPoint((A.x() + B.x()) / 2.0,  (A.y() + B.y()) / 2.0, (A.z() + B.z()) / 2.0);

	CGALSphere smallestCircumsphere(sphereCenter, pow(sphereRadius, 2));
	
	float encroachingCandidateDistance;
	vector<float> circumradiusMap; // value is computed only for those points which are encroaching
	
	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
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

	for (LCC::One_dart_per_cell_range<0>::iterator pointIter = plc.one_dart_per_cell<0>().begin(), pointIterEnd = plc.one_dart_per_cell<0>().end(); pointIter != pointIterEnd; pointIter++, vertexId++)
		if (circumradiusMap[vertexId] != INVALID_VALUE)
			if (maxCircumradius < circumradiusMap[vertexId])
			{
				maxCircumradius = circumradiusMap[vertexId];
				refPoint = &(plc.point(pointIter));
			}	

	return;
}


/*! \fn float dotProduct(DartHandle segment1Handle, DartHandle segment2Handle)
    \brief Computes dot product of vectors represented by input segments.

    \param [in] segment1Handle DartHandle for the first segment
    \param [in] segment2Handle DartHandle for the second segment
 */
float dotProduct(DartHandle segment1Handle, DartHandle segment2Handle)
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

/*! \fn float vectorMagnitude(DartHandle inputSegmentHandle)
    \brief Computes magnitue of input vector

    \param [in] inputSegmentHandle DartHandle of input segment
 */
float vectorMagnitude(DartHandle inputSegmentHandle)
{
	CGALPoint segmentVertices[2];

	segmentVertices[0] = plc.point(inputSegmentHandle);
	segmentVertices[1] = plc.point(plc.beta(inputSegmentHandle, 1));	

	CGALPoint vector = CGALPoint(segmentVertices[0].x() - segmentVertices[1].x(), segmentVertices[0].y() - segmentVertices[1].y(), segmentVertices[0].z() - segmentVertices[1].z());

	float vectorMagnitude = sqrtf(powf(vector.x(), 2.0) + powf(vector.y(), 2.0) + powf(vector.z(), 2.0));

	return vectorMagnitude;
}


/*! \fn float computeAngleBetweenSegments(DartHandle segment1Handle, DartHandle segment2Handle)
    \brief Computes angle between vectors represented by input segments.

    \param [in] segment1Handle DartHandle for first input segment
    \param [in] segment2Handle DartHandle for second input segment
 */
float computeAngleBetweenSegments(DartHandle segment1Handle, DartHandle segment2Handle)
{

	float angle = acosf(dotProduct(segment1Handle, segment2Handle) / (vectorMagnitude(segment1Handle) * vectorMagnitude(segment2Handle)));

	return (angle * 180.0f / Pi); // conversion to degrees
}

/*! \fn bool isVertexAcute(DartHandle inputPointHandle)
    \brief Tests whether input vertex is acute

    A vertex is _not_ acute if there exists _atleast_ one pair of _incident_ segments which have angle _greater_ than 90 degrees.    
   \param [in] inputPointHandle DartHandle of input vertex
 */
bool isVertexAcute(DartHandle inputPointHandle)
{
	// Determine segment-pair(involving A) 
	vector<DartHandle> incidentOnInputPoint;

	 
	for (LCC::One_dart_per_incident_cell_range<1, 0>::iterator incidentSegmentIter = plc.one_dart_per_incident_cell<1, 0>(inputPointHandle).begin(), incidentSegmentIterEnd = plc.one_dart_per_incident_cell<1, 0>(inputPointHandle).end(); incidentSegmentIter != incidentSegmentIterEnd;incidentSegmentIter++)
		incidentOnInputPoint.push_back(incidentSegmentIter);

	// Compute angle between all possible pairs(NAIVE SOLUTION)
	for (vector<DartHandle>::iterator segIter1 = incidentOnInputPoint.begin(); segIter1 != incidentOnInputPoint.end(); segIter1++) 
		for (vector<DartHandle>::iterator segIter2 = incidentOnInputPoint.begin(); *segIter1 != *segIter2 && segIter2 != incidentOnInputPoint.end(); segIter2++)
			if (computeAngleBetweenSegments(*segIter1, *segIter2) < 90.0f)
				return true;

	return false; 
}


/*! \fn unsigned int determineSegmentType(DartHandle missingSegmentHandle)
    \brief Determines the input segment type

    A segment is of _type-1_ if _none_ of its endpoints are _acute_. If _exactly_ one endpoint is acute than segment is of _type-2_.
    \param [in] DartHandle for input constraint segment

 */
unsigned int determineSegmentType(DartHandle missingSegmentHandle)
{

	CGALPoint &A = plc.point(missingSegmentHandle);
	CGALPoint &B = plc.point(plc.beta(missingSegmentHandle, 1));

	bool vertexAIsAcute = isVertexAcute(missingSegmentHandle);
	bool vertexBIsAcute = isVertexAcute(plc.beta(missingSegmentHandle, 1));

	if (!vertexAIsAcute && !vertexBIsAcute)	
		return 1;
	
	else if (vertexAIsAcute != vertexBIsAcute) // effectively XOR
		return 2;
	else
		return 3;
}

/*! \fn float computeSegmentLangth(CGALPoint &A, CGALPoint &B)
    \brief Computes Euler length of the segment represented by _A_ and _B_

    \param [in] A First endpoint of the segment
    \param [in] B Second endpoint of the segment
 */
float computeSegmentLength(CGALPoint &A, CGALPoint &B)
{
	float sLength;
	sLength = sqrt(pow(A.x() - B.x(), 2) + pow(A.y() - B.y(), 2) + pow(A.z() - B.z(), 2));
	return sLength;
}

/*! void updatePLCAndDT(CGALPoint &v, DartHandle missingSegmentHandle)
   /brief Updates PLC and Delaunay triangulation after insertion of a new vertex into a constraint segment.
    
   /param [in] v New vertex to be inserted into a constraint segment of PLC
   /param p[ih] missignSegmentHandle DartHandle of constraint segment which will be split after inserting new vertex  
 */
void updatePLCAndDT(CGALPoint &v, DartHandle missingSegmentHandle)
{

	// update PLC
	plc.insert_point_in_cell<1>(missingSegmentHandle, v);
	// update DT
	computeDelaunayTetrahedralization(); 
}

/*! \fn void splitMissingSegment(DartHandle missingSegmentHandle)
    \brief Splits the missing constraint segment
    
    A constraint segment is split(by insertion of a new vertex) by applying segment splitting rules. Purpose of splitting is to make resulting subsegments _strongly Delaunay_.
    \param [in] missingSegmentHandle DartHandle of missing constraint segment
 */
void splitMissingSegment(DartHandle missingSegmentHandle)
{

	CGALPoint vb, refPoint;
	CGALPoint sphereCenter;
	float sphereRadius;
	unsigned int segmentType;
	segmentType = determineSegmentType(missingSegmentHandle);

	computeReferencePoint(&refPoint, missingSegmentHandle);
	
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

		CGALSphericalSegment seg(p1, p2);
		CGALSphericalLineArc lineArc(seg);
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
			CGALPoint t;
			t = A;
			A = B;
			B = t;
		}
		
		
		// Segment calculations
		
		ApB[0] = acuteParentHandle;
		ApB[1] = BHandle;
		
	
		CGALSphericalPoint p1(plc.point(acuteParentHandle).x(), plc.point(acuteParentHandle).y(), plc.point(acuteParentHandle).z());
		CGALSphericalPoint p2(plc.point(BHandle).x(), plc.point(BHandle).y(), plc.point(BHandle).z());
		
		CGALSphericalSegment seg(p1, p2);
		CGALSphericalLineArc lineArc(seg);
		
		/// Sphere calculations
		CGALSphericalPoint acuteParent(plc.point(acuteParentHandle).x(), plc.point(acuteParentHandle).y(), plc.point(acuteParentHandle).z());
		
		CGALPoint acuteParentLinearField(plc.point(acuteParentHandle).x(), plc.point(acuteParentHandle).y(), plc.point(acuteParentHandle).z());
	
		float ApRefPointLength = computeSegmentLength(acuteParentLinearField, refPoint);

		CGALSphericalSphere s(acuteParent, pow(ApRefPointLength, 2));

		vector<Object> intersections;

		CGAL::intersection(lineArc, s, back_inserter(intersections));
	
		if (intersections.size() > 0)
		{	
			v = object_cast<CGALPoint>(intersections.back());
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
			CGALSphericalPoint sphereCenter(acuteParent);
			unsigned int avLength = computeSegmentLength(plc.point(AHandle), v);
			if (vrefpointLength < 0.5 * avLength)
				{
					CGALPoint temp(acuteParentLinearField.x(), acuteParentLinearField.y(), acuteParentLinearField.z());
					acuteparentALength = computeSegmentLength(temp, plc.point(AHandle));
					sphereRadius = acuteparentALength + avLength - vrefpointLength;	
				}
			else
				sphereRadius = acuteparentALength + 0.5 * avLength;
		        
			s = CGALSphericalSphere(sphereCenter, pow(sphereRadius, 2));
		
			CGAL::intersection(lineArc, s, back_inserter(intersections));

			if (intersections.size() > 0)
			{
				v = object_cast<CGALPoint>(intersections.back());
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

/*! void recoverConstraintSegments()
    \brief Top-level routine for recovering constraint segments.
	
    Missing constriant segments of PLC are recovered by spliting the segments in such a way so that resulting sub-segments become _strongly Delaunay_.
*/
void recoverConstraintSegments()
{
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

/* \class DegenerateVertexSetCandidate
   \brief Models a set of 5 _cospherical_ vertices   
 */
class DegenerateVertexSetCandidate
{
	public:
		DartHandle pointHandles[5]; /*< Array of DartHandles of _cospherical_ vertices*/
};


/*! bool areCospherical(DegenerateVertexSetCandidate degenSet)
    \brief Tests whether input 5-vertex set is _cospherical_

    \param [in] degenSet A set of 5 vertices
*/
bool areCospherical(DegenerateVertexSetCandidate degenSet)
{	
	CGALPoint p[5];
	
	for (unsigned int i = 0; i < 5; i++)	
		p[i] = plc.point(degenSet.pointHandles[i]);

	if (side_of_bounded_sphere(p[0], p[1], p[2], p[3], p[4]) == ON_BOUNDARY)
		return true;
	else
		return false;
}


/*! \fn void addLocalDegeneraciesToQueue(vector<DegenerateVertexSetCandiate> &localDegeneracySet)
    \brief Collects all local degneracies from PLC into a queue.
    
     \param [out] Write local degeneracy sets in this queue    
 */
void addLocalDegeneraciesToQueue(vector<DegenerateVertexSetCandidate> &localDegeneracySet)
{	
	DegenerateVertexSetCandidate degenerateSetCandidate;
	
	for (Delaunay::Finite_cells_iterator cellIter = DT.finite_cells_begin(); cellIter != DT.finite_cells_end(); cellIter++)
	{
		for (unsigned int n = 0; n < 3; n++)
			degenerateSetCandidate.pointHandles[n] = (cellIter->vertex(n))->info(); // info structure contains pointIds

		for (unsigned int j = 0; j < 4; j++)
			{
				VertexHandle vh = DT.mirror_vertex(cellIter, j);	
				degenerateSetCandidate.pointHandles[4] = vh->info();		
						
				if (areCospherical(degenerateSetCandidate))
				{
					localDegeneracySet.push_back(degenerateSetCandidate);	
				}
			}
	}
}

			
/* \fn bool isVertexPerturbable(unsigned int pointId)
   \brief Tests whether input vertex is _perturbable_

   \param [in] pointHandle DartHandle for the input point
 */
bool isVertexPerturbable(DartHandle pointHandle)
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

/*! \fn bool isVertexSegmentSafePerturbable(DartHandle pointHandle)

    \brief Tests whether input vertex is _segment safe_ perturbable  	
    
    A vertex is _segment-safe_ perturbable if there exists a perturbation of vertex which does not make PLC _inconsistent_.
    
    \param [in] pointHandle DartHandle of input point
*/
bool isVertexSegmentSafePertubable(DartHandle pointHandle)
{
	bool segmentSafePerturbable = false;

	return segmentSafePerturbable;
}

/*! \fn bool isDegeneracyRemovable(DegeneracteVertexSetCandidate degenCandidate)
    \brief Tests whether input degeneracy is _removable_
    
    A degeneracy is _removable_ if it has atleast _one_ vertex with a possibility of _segment-safe_ perburbation.    

    \param [in] degenCandidate A degenerate vertex set candidate to be tested for _removability_.
 */
bool isDegeneracyRemovable(DegenerateVertexSetCandidate degenCandidate)
{
	bool removable = false;
	return removable;
}

/* /fn void perturbRemove(DartHandle pointHandle, vector<DegenerateVertexSetCandidate>& localDegeneracySet)
   /brief Removes vertex using _perturbation_ for breaking local degeneracy.   
 
*/
void perturbRemove(unsigned int setIndex, vector<DegenerateVertexSetCandidate>& localDegeneracySet)
{
	// DO NOTHING FOR NOW
}

/*! \fn bool areAffinelyIndependent(DegenerateVertexSetCandidate degenSet)

    \brief Tests whether vertices in the input local deheneracy set have _affine_ independence.

    \param [in] degenSet A set of 5 _cospherical_ vertices
 */
bool areAffinelyIndependent(DegenerateVertexSetCandidate degenSet)
{

	CGALPoint v1 = plc.point(degenSet.pointHandles[0]);
	CGALPoint v2 = plc.point(degenSet.pointHandles[1]);
	CGALPoint v3 = plc.point(degenSet.pointHandles[2]);
	CGALPoint v4 = plc.point(degenSet.pointHandles[3]);
	CGALPoint v5 = plc.point(degenSet.pointHandles[4]);

	if (coplanar(v1, v2, v3, v4) || coplanar(v2, v3, v4, v5) || coplanar(v3, v4, v5, v1) || coplanar(v4, v5, v1, v2) || coplanar(v5, v1, v2, v3))
		return false;

	return true; 
}


/*! \fn void computeBreakPoint(CGALPoint &vb, unsigned int pointId, vector<DegenerateVertexSetCandidate> &localDegeneracySet)
    
    \brief Computes _breaking point_ used for which is used for resolving local degeneracy.

    \param vb [out] The computed breaking point

    \param pointId [in] Index into the input local degeneracy set

    \param localDegeneracySet [in] The local degeneracy set containing all local degeneracies of PLC
 */
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


/*! \fn bool isVertexEncroachingSegment(CGALPoint vb, DartHandle segmentHandle)
 
    \brief Tests whether input vertex _encroaches_ the segment.
    
    A vertex _encroaches_ a segment if it lies _on_ or _inside_ its diameter sphere.

    \param [in] vb Input vertex

    \param [in] DartHandle to the segment to be tested for _encroachment_
*/
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

/*! \fn bool isVertexEncroachingFace(CGALPoint vb, DartHandle facetHandle)
 
    \brief Tests whether input vertex _encroaches_ the facet.
    
    A vertex _encroaches_ a facets if it lies _on_ or _inside_ its diameter circumsphere.

    \param [in] vb Input vertex

    \param [in] DartHandle to the facet to be tested for _encroachment_
*/
bool isVertexEncroachingFace(CGALPoint vb, DartHandle facetHandle)
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


/*! \fn bool isEncroachingPLC(CGALPoint vb)
    
    \brief Tests whether input point encroaches any segment or a facet of PLC.  
    
    \param [in] vb Input point
 */
bool isEncroachingPLC(CGALPoint vb)
{
	// encroaches any segment
	for (LCC::One_dart_per_cell_range<1>::iterator segmentIter = plc.one_dart_per_cell<1>().begin(), segmentIterEnd = plc.one_dart_per_cell<1>().end(); segmentIter != segmentIterEnd; segmentIter++)
		if (isVertexEncroachingSegment(vb, segmentIter))
			return true;

	// encroaches any face
	for (LCC::One_dart_per_cell_range<2>::iterator faceIter = plc.one_dart_per_cell<2>().begin(), faceIterEnd = plc.one_dart_per_cell<2>().end(); faceIter != faceIterEnd; faceIter++)
		if (isVertexEncroachingFace(vb, faceIter))
			return true;

	return false;
}


/*! \fn void boundaryProtection(CGALPoint vb)

  \brief Boundary protection procedure

  \param [in] vb Input point
 */
void boundaryProtection(CGALPoint vb)
{

	for (LCC::One_dart_per_cell_range<1>::iterator segmentIter = plc.one_dart_per_cell<1>().begin(), segmentIterEnd = plc.one_dart_per_cell<1>().end(); segmentIter != segmentIterEnd; segmentIter++)
		if (isVertexEncroachingSegment(vb, segmentIter))
		{
			// circumcenter of segment 'n':
			CGALPoint endpoint1(plc.point(segmentIter));
		        CGALPoint endpoint2(plc.point(plc.beta(segmentIter, 1)));	
			float x = (endpoint1.x() + endpoint2.x()) / 2.0;
			float y = (endpoint1.y() + endpoint2.y()) / 2.0;
			float z = (endpoint1.z() + endpoint2.z()) / 2.0;

			CGALPoint sphereCenter(x, y, z);

			updatePLCAndDT(sphereCenter, segmentIter);
		}
	
}


/*! \fn void removeLocalDegeneracies()
    \brief Top-level function for dealing with removal of local degeneracies from PLC.
 */
void removeLocalDegeneracies()
{

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
					plc.create_dart(vb); // TODO: Is it connected with rest of the cell complex?
			}						
		}
	}
	cout << "\nLocal degeneracy removal completed";
}

/*! \fn void formMissingSubfaceQueue(vector<DartHandle> &missingSubfacesQueue)
    
     \brief Collects all constraint faces from PLC which are missing from corresponding Delaunay triangulation.

     \param [out] missingSubfacesQueue Add all missing constraint faces in this vector
 */
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

/*! \fn void formCavity(vector<DartHandle> *cavity, DartHandle missingSubfaceHandle, vector<DartHandle>& lcc3CellsToBeRemoved)

     \brief Computes cavities created by insertion of missing constraint faces into Delaunay triangulation.

     \param [out] cavity Output cavity

     \param [in] DartHandle for a missing constraint facet

     \param [out] Collect handles to the 3-cells intersecting input constraint facet
 */
void formCavity(vector<DartHandle> *cavity, DartHandle missingSubfaceHandle, vector<DartHandle>& lcc3CellsToBeRemoved)
{
	vector<DartHandle> intersectingTets;				
	CGALPoint pTet[4], pTri[3];
	
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
	for (LCC::One_dart_per_cell_range<2>::iterator faceIter = cdtMesh.one_dart_per_cell<2>().begin(), faceIterEnd = cdtMesh.one_dart_per_cell<2>().end(); faceIter != faceIterEnd; faceIter++)
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
	
	CGALPoint v1 = plc.point(missingSubfaceHandle); 
	CGALPoint v2 = plc.point(plc.beta(missingSubfaceHandle, 1));
	CGALPoint v3 = plc.point(plc.beta(missingSubfaceHandle, 1, 1));	

	// generating random points
	Random rnd1(30.0), rnd2(10.0);

	float s = rnd1.uniform_real<float>(); 
	float t = rnd2.uniform_real<float>();
	CGALPoint points[3];
	
	float xRand, yRand, zRand;

	for (unsigned int g = 0, n = 0; g < globalCavity.size(); g++)
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


/*! \fn bool isFaceStronglyDelaunay(DartHandle facetHandle, vector<DartHandle> cavityVerticesSet, LCC cavityLCC)
    
    \brief Tests whether the input facet is _strongly Delaunay_ with respect to cavity vertices.

    \param [in] facetHandle The input facet

    \param [in] cavityVerticesSet Set of vertices in the cavity 

    \param [in] cavityLCC Linear cell complex representation of cavity
 */
bool isFaceStronglyDelaunay(DartHandle facetHandle, vector<DartHandle> cavityVerticesSet, LCC cavityLCC)
{
	// compute Delunay tetrahedralization of vertices
		// check if the face is there in that DT 
		// If yes then check if the enclosing sphere has any other vertex on its surface
			// If yes, then it is Delaunay but not strongly Delaunay, return false
			// If no, then it is strongly Delaunay, return true
		// If no, then it is not even Delaunay, return false 	
		
	Delaunay tempDT;
	vector<CGALPoint> cavityVertices;
	bool faceStronglyDelaunay;
	
	for (vector<DartHandle>::iterator vertexIter = cavityVerticesSet.begin(), vertexIterEnd = cavityVerticesSet.end(); vertexIter != vertexIterEnd; vertexIter++)
		cavityVertices.push_back(cavityLCC.point(*vertexIter));

	// Compute DT
	tempDT.insert(cavityVertices.begin(), cavityVertices.end());

	// TODO: Here we only need to check if input facet is present in DT becuase local degeneracies are implicitly removed in CGAL's DT algorithm
	VertexHandle vh1, vh2, vh3;
	CellHandle ch;
	int i, j, k;

	for (LCC::One_dart_per_cell_range<2>::iterator facetIter = plc.one_dart_per_cell<2>().begin(), facetIterEnd = plc.one_dart_per_cell<2>().end(); facetIter != facetIterEnd; facetIter++)
	{
		tempDT.is_vertex(plc.point(facetIter), vh1); 
		tempDT.is_vertex(plc.point(plc.beta(facetIter, 1)), vh2);
		tempDT.is_vertex(plc.point(plc.beta(facetIter, 1, 1)), vh3);
		if (tempDT.is_facet(vh1, vh2, vh3, ch, i, j, k))
			faceStronglyDelaunay = true;
		else
			faceStronglyDelaunay = false;
	}
	
	return faceStronglyDelaunay;
}


/*! \fn bool areSameFacets(LCC lcc1, DartHandle d1, LCC lcc2, DartHandle d2)

    \brief Tests whether two input facets from different LCC's are same.

    \param [in] lcc1 LCC containing first facet

    \param [in] d1 DartHandle for first facet

    \param [in] lcc2 LCC constaining second facet

    \param [in] d2 DartHandle for second facet
 */
bool areSameFacets(LCC lcc1, DartHandle d1, LCC lcc2, DartHandle d2)
{
	DartHandle d3 = lcc2.make_triangle(lcc1.point(d1), lcc1.point(lcc1.beta(d1, 1)), lcc1.point(lcc1.beta(d1, 1, 1))); // TODO: Fix this roundabout way of comparing facets from 2 LCCs
	bool areSameFacets = lcc2.are_facets_same_geometry(d3, d2);
	remove_cell<LCC, 2>(lcc2, d3);	

	return areSameFacets;
}

/*! \fn DartHandle locateFacetInCavity(DartHandle facet, LCC cavityLCC)

    \brief Returns dart handle of the facet from cavityLCC corresponding to a facet in cdtMesh.

    \param [in] cdtMeshFacetHandle DartHandle to a facet in cdtMesh

    \param [in] cavityLCC LCC representation of the cavity
 */
DartHandle locateFacetInCavity(DartHandle cdtMeshFacetHandle, LCC cavityLCC)
{
	for (LCC::One_dart_per_cell_range<2>::iterator facetIter = cavityLCC.one_dart_per_cell<2>().begin(), facetIterEnd = cavityLCC.one_dart_per_cell<2>().end(); facetIter != facetIterEnd; facetIter++)
	{
		if (areSameFacets(cdtMesh, cdtMeshFacetHandle, cavityLCC, facetIter))
			return facetIter;
		else
			continue;
	}
}


/*! \fn void addFaceToLCC(CGALPoint p1, CGALPoint p2, CGALPoint p3, LCC &cavityLCC)

    \brief Adds a facet defined by p1, p2, p3 to the cavity.

    \param [in] p1 First point

    \param [in] p2 Second point

    \param [in] p3 Third point

    \param [out] cavityLCC Add facet defined by p1, p2, p3 in this cavity
 */
void addFaceToLCC(CGALPoint p1, CGALPoint p2, CGALPoint p3, LCC &cavityLCC)
{
	DartHandle newFace = cavityLCC.make_triangle(p1, p2, p3);
	int sewedMark = cavityLCC.get_new_mark();
	vector<pair<DartHandle, DartHandle> > twoCellsToBeSewed;
	
	// sew the triangle with existing faces in cavityLCC
	for (LCC::One_dart_per_incident_cell_range<1, 2>::iterator segIter1 = cavityLCC.one_dart_per_incident_cell<1, 2>(newFace).begin(), segIterEnd1 = cavityLCC.one_dart_per_incident_cell<1, 2>(newFace).end(); segIter1 != segIterEnd1; segIter1++)
	{
		if (!cavityLCC.is_marked(segIter1, sewedMark)) // not sewed till now
		{
			for (LCC::Dart_range::iterator segIter2 = plc.darts().begin(), segIterEnd2 = plc.darts().end(); segIter2 != segIterEnd2; segIter2++)
			{
				if (!plc.is_marked(segIter2, sewedMark) && plc.is_sewable<2>(segIter1, segIter2))
				{
					if (areGeometricallySameSegments(segIter1, segIter2)) // checks the geometry of segments
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

}


namespace std
{
	template<> struct hash<CGALPoint>
	{
		size_t operator () (const CGALPoint& p) const
		{
			return sqrt(pow(p.x(), 2) + pow(p.y(), 2) + pow(p.z(), 2)); // TODO: effective hash
		}
	};
}


/*! \fn bool isTetOutsideCavity(DartHandle cellHandle, LCC filledCavityLCC)

    \brief Tests whether a tetrahedron is outside cavity using Point in polygon approach.
   
    \param [in] cellHandle DartHandle to a 3-cell 
    
    \param [in] filledCavityLCC LCC representation of the tetrahedralized cavity
 */
bool isTetOutsideCavity(DartHandle cellHandle, LCC filledCavityLCC)
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
				
		CGALPoint p[4];

		unsigned int i = 0;
		for (LCC::One_dart_per_incident_cell_range<0, 3>::iterator vertexIter = filledCavityLCC.one_dart_per_incident_cell<0, 3>(cellHandle).begin(), vertexIterEnd = filledCavityLCC.one_dart_per_incident_cell<0, 3>(cellHandle).end(); vertexIter != vertexIterEnd; vertexIter++)
			p[i++] = filledCavityLCC.point(vertexIter);
		
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

		CGALPoint tetBarycenter(x, y, z);
	
		// generate a random ray
		
		Random rndGenerator1(50.0); // some numbers given as seed
		Random rndGenerator2(60.0);
		Random rndGenerator3(70.0);

		CGALPoint randomPoint(rndGenerator1.uniform_real<float>(), rndGenerator2.uniform_real<float>(), rndGenerator3.uniform_real<float>());

		CGALRay randomRay(tetBarycenter, randomPoint); // ray starting at tetBarycenter and passing through randomPoint
		unsigned int intersectionCount = 0;
		
		bool isUniqueIntersectionPoint = false;
		unordered_set<CGALPoint> uniqueIntersectionPoints;
		for (LCC::One_dart_per_cell_range<2>::iterator facetIter = filledCavityLCC.one_dart_per_cell<2>().begin(), facetIterEnd = filledCavityLCC.one_dart_per_cell<2>().end(); facetIter != facetIterEnd; facetIter++)
		{
			// compute its intersection with all faces of cavity
			// count number of intersections
			CGALPoint facetVertices[3];
			
			unsigned int k = 0;	
			for (LCC::One_dart_per_incident_cell_range<0, 2>::iterator pointIter = filledCavityLCC.one_dart_per_incident_cell<0, 2>(facetIter).begin(), pointIterEnd = filledCavityLCC.one_dart_per_incident_cell<0, 2>(facetIter).end(); pointIter != pointIterEnd; pointIter++)
				facetVertices[k++] = filledCavityLCC.point(pointIter);
		

			CGALTriangle cavityFacet(facetVertices[0], facetVertices[1], facetVertices[2]);
		        cpp11::result_of<K::Intersect_3(CGALRay, CGALTriangle)>::type intersectionResult;
			CGALPoint *intersectionPoint;

			if ((intersectionResult = intersection(randomRay, cavityFacet)) != NULL)
			{
				intersectionPoint = boost::get<CGALPoint>(&*intersectionResult);
				isUniqueIntersectionPoint = (uniqueIntersectionPoints.insert(*intersectionPoint)).second;	
				if (isUniqueIntersectionPoint)
					intersectionCount++;
				else
					continue;
			}
			else
				continue;
		}
		
		if (intersectionCount % 2 == 0)
			return true;
		else
			return false;	
}               

/*! \fn void createLCCRepresentationForCavity()
    \brief Constructs & returns LCC representation for the input cavity
    \param [in] cavity Cavity represented as vector of DartHandles over cdtMesh
    \param [out] cavityLCC LCC representation of cavity, contains DartHandle to corresponding face in cdtMesh as its info<2>
    \param [out] Mapping between 2-cells in cavityLCC and cdtMesh    
*/
void createLCCRepresentationForCavity(vector<DartHandle> cavity, LCCWithInfo &cavityLCC)
{
	for (vector<DartHandle>::iterator faceIter = cavity.begin(), faceIterEnd = cavity.end(); faceIter != faceIterEnd; faceIter++)
	{
		CGALPoint p[3];
		p[0] = plc.point(*faceIter);
		p[1] = plc.point(plc.beta(*faceIter, 1));
		p[2] = plc.point(plc.beta(*faceIter, 1, 1));

		DartHandleWithInfo d = cavityLCC.make_triangle(p[0], p[1], p[2]);
		for(LCCWithInfo::Dart_of_cell_range<2, 3>::iterator dartIter = cavityLCC.dart_of_cell<2, 3>(d), dartIterEnd = cavityLCC.dart_of_cell<2, 3>(d).end(); dartIter != dartIterEnd; vIter++)
			cavityLCC.info<0>(dartIter) = *faceIter; // TODO: Correct??
	}
		// TODO: sew the face 
}


/*! \fn void cavityRetetrahedralization(vector <DartHandle>& cavity, vector<DartHandle>& lcc3CellsToBeRemoved)

  \brief Retetrahedralizes the cavity created due to removal of 3-cells intersecting the missing constraint face.

  First it expands the cavity further(if required) to ensure all facets of cavity are _strongly Delaunay_. Then it fills the cavity with 3-cells and sews it with the original mesh to generate mesh containing the missing constraint face.

  \param [out] cavity Vector of DartHandles for faces representing cavity

  \param [out] Vector of DartHandles for 3-cells to be removed 
 */
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
	LCCWithInfo cavityLCC;
	createLCCRepresentationForCavity(cavity, cavityLCC);  
	
	
	vector<DartHandle> nonStronglyDelaunayFacesInCavity;
	vector<DartHandle> cavityVerticesSet; // V

	// compute set of non-strongly Delunay faces in cavity
	do
	{	
		for (LCCWithInfo::One_dart_per_cell_range<0>::iterator vertexIter = cavityLCC.one_dart_per_cell<0>().begin(), vertexIterEnd = cavityLCC.one_dart_per_cell<0>().end(); vertextIter != vertexIterEnd; vertexIter++)
			cavityVerticesSet.push_back(vertexIter);

		for (LCCWithInfo::One_dart_per_cell_range<2>::iterator facetIter = cavityLCC.one_dart_per_cell<2>().begin(), facetIterEnd = cavityLCC.one_dart_per_cell<2>().end(); facetIter != facetIterEnd; facetIter++)
		{
			if (isFaceStronglyDelaunay(facetIter, cavityVerticesSet, cavityLCC))
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
			DartHandleWithInfo cellToBeRemovedHandle = cavityLCC.info<2>(tempNonStronglyDelaunayFacet);
			lcc3CellsToBeRemoved.push_back(cellToBeRemovedHandle);

			for (LCC::One_dart_per_incident_cell_range<2, 3>::iterator faceIter = cdtMesh.one_dart_per_incident_cell<2, 3>(cellToBeRemovedHandle).begin(), faceIterEnd = cdtMesh.one_dart_per_incident_cell<2, 3>(cellToBeRemovedHandle).end(); facetIter != facetIterEnd; faceIter++) // for all faces of the removed cell
			{
				if (!areSameFacet(faceIter, tempNonStronglyDelaunayFacet))
				{
					if ((facetHandle = locateFacetInCavity(faceIter, cavityLCC)) != -1)
						cavityLCC.remove_cell<2>(facetHandle);
					else // add facet to the cavity
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
			
	// create a vector of vertices of cavity
	vector <CGALPoint> cavityVertices; 
	for (LCCWithInfo::One_dart_per_cell<0>::iterator vertexIter = cavityLCC.one_dart_per_cell<0>().begin(), vertexIterEnd = cavityLCC.one_dart_per_cell<0>().end(); vertexIter != vertexIterEnd; vertexIter++)
		cavityVertices.push_back(cavityLCC.point(vertexIter));

	Delaunay cavityDT;

	cavityDT.insert(cavityVertices.begin(), cavityVertices.end()); 

	LCC filledCavityLCC;

	import_from_triangulation_3(cavityDT, filledCavityLCC);	
	
	// Mark tetrahedrons as 'inside/outside' wrt. cavity boundary	
	
	// Remove outside tetrahedrons
	// Sew the remaining 3-cells 
	int insideOutsideMark = filledCavityLCC.get_new_mark();

	for (LCC::One_dart_per_cell_range<3>::iterator tetIter = filledCavityLCC.one_dart_per_cell<3>().begin(), tetIterEnd = filledCavityLCC.one_dart_per_cell<3>().end(); tetIter != tetIterEnd; tetIter++)
	{
		if (isTetOutsideCavity(tetIter, filledCavityLCC))
			filledCavityLCC.mark(tetIter, insideOutsideMark);
		else
			continue;
	}
	
	// remove the marked 3-cells & add remaining ones to cdtMesh
	for (LCC::One_dart_per_cell_range<3>::iterator removeTetIter = filledCavityLCC.one_dart_per_cell<3>().begin(), removeTetIterEnd = filledCavityLCC.one_dart_per_cell<3>().end(); removeTetIter != removeTetIterEnd; removeTetIter++)
	{
		if(filledCavityLCC.is_marked(removeTetIter))
			filledCavityLCC.remove_cell<3>(removeTetIter);
		else
		{
			CGALPoint p[4];
			unsigned int i = 0;
			for (LCC::One_dart_per_incident_cell_range<0, 3>::iterator vIter = filledCavityLCC.one_dart_per_incident_cell<0, 3>(removeTetIter).begin(), vIterEnd = filledCavityLCC.one_dart_per_incident_cell<0, 3>(removeTetIter).end(); vIter != vIterEnd; vIter++)
				p[i++] = filledCavityLCC.point(vIter);
			
			cdtMesh.make_tetrahedron(p[0], p[1], p[2], p[3]);
		}
	}
	
	// sew the new cavity filled LCC with existing cdtMesh
	cdtMesh.sew3_same_facets();
		
}


/*! \fn void recoverConstraintFaces()

    \brief Top-level routine for recovering constraint faces.
 */
void recoverConstraintFaces()
{
	 			
        vector<DartHandle> missingSubfacesQueue;
	formMissingSubfaceQueue(missingSubfacesQueue);
	vector<DartHandle> cavity[2];
	unsigned int missingSubfaceId;

	vector<Triangle> cdtFacetList;
	vector<DartHandle> lcc3CellsToBeRemoved;

	import_from_triangulation_3(cdtMesh, DT);

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
	

/*! Main procedure
 */
int main()
{
	readPLCInput();
	computeDelaunayTetrahedralization();
	recoverConstraintSegments();
	removeLocalDegeneracies();
	recoverConstraintFaces();

	return 0;
}

