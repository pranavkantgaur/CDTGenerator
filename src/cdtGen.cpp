#include "cdtGen.h"

static size_t pointId = 0;
static float tempPoint[3];
static size_t dimensionId = 0;
vector<CGALPoint> plcVertexVector;
vector<Triangle> plcFaceVector;

/*! \fn static int vertex_cb(p_ply_argument argument)
    \brief Callback for reading vertex from PLY file	
    \param [in] argument represents PLY file
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
    \param [in] argument Represents PLY file 
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

/*! \fn bool CDTGenerator::areGeometricallySameSegments(DartHandle d1, DartHandle d2)
    \brief Tests whether segments represented by d1 and d2 are _geometrically_ same.
    \param [in] d1 DartHandle for first segment
    \param [in] d2 DartHandle for second segment
*/
bool CDTGenerator::areGeometricallySameSegments(DartHandle d1, DartHandle d2)
{
	if (plc.point(d1) == plc.point(plc.beta(d2, 1)))
		if (plc.point(plc.beta(d1, 1)) == plc.point(d2))
			return true;
	return false;
} 


/*! \fn void CDTGenerator::readPLCInput()
    \brief Reads PLC from input(only PLY supported currently) file.

    Reads vertex and face information and initializes corresponding linear cell complex _plc_.
*/
void CDTGenerator::readPLCInput()
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
	dimensionId = 0;

	if (!ply_read(inputPLY))
	{
		cout << "Cannot read the PLY file :(";
		exit(0);
	}

	ply_close(inputPLY);

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
}


/*! void CDTGenerator::markInfinteVertexDart(LCC::Dart_handle d)
*   \brief Marks all darts defining the infinite vertex in LCC.
*   \param [in] d Handle to the dart representing infinite vertex.
*   \param [in] lcc Linear cell complex to be marked.
*   \param [in] infiniteVertexMark Mark representing a dart defining the infinite vertex.
*/
void CDTGenerator::markInfiniteVertexDart(LCCWithIntInfo::Dart_handle d, LCCWithIntInfo &lcc, int infiniteVertexMark)
{
	if (!lcc.is_marked(d, infiniteVertexMark))
	{
		lcc.mark(d, infiniteVertexMark);
		
		d = lcc.beta(d, 2, 1);
		markInfiniteVertexDart(d, lcc, infiniteVertexMark);

		d = lcc.beta(d, 3, 1);
		markInfiniteVertexDart(d, lcc, infiniteVertexMark);
	}

	return;
}


/** \fn isInfinite(LCC::Dart_handle adart, LCC lcc, unsigned int cell_dimension)
 *  \brief Tests whether the given i-cell is infinite.
 *  \param [in] adart a dart handle to the i-cell.
 *  \param [in] lcc Linear cell complex to be checked.
 *  \param [in] infiniteVertexMark mark representing the infinite vertex.
 *  \param [in] cell_dimension dimension of i-cell.
 *  \return True if input vertex or face is infinite    
**/
bool CDTGenerator::isInfinite(LCCWithIntInfo::Dart_handle adart, LCCWithIntInfo lcc, int infiniteVertexMark, unsigned int cell_dimension)
{
	bool isInfinite = false;
	
	if (cell_dimension == 0)
	{
		if (lcc.is_marked(adart, infiniteVertexMark)) 
		return true;
	}

	if (cell_dimension == 2)
	{
		for (LCCWithIntInfo::One_dart_per_incident_cell_range<0, 2>::iterator pIter = lcc.one_dart_per_incident_cell<0, 2>(adart).begin(), pIterEnd = lcc.one_dart_per_incident_cell<0, 2>(adart).end(); pIter != pIterEnd; pIter++)
		{
			if (lcc.is_marked(pIter, infiniteVertexMark)) 
			{
				isInfinite = true;
				break;
			}
		}
	}
	return isInfinite;
}


/*! \fn void CDTGenerator::writePLYOutput(LCC::Dart_handle dartToInfiniteVertex, LCC &lcc, string fileName)
*   \brief Writes mesh represented as LCC to PLY file
*   \param [in] dartToInfiniteVertex Dart handle to the infinite vertex.
*   \param [in] lcc Linear cell complex.
*   \param [in] fileName Name of output PLY file.
*/
void CDTGenerator::writePLYOutput(LCCWithIntInfo::Dart_handle dartToInfiniteVertex, LCCWithIntInfo &lcc, string fileName)
{
	p_ply lccOutputPLY;

	if ((lccOutputPLY = ply_create(fileName.c_str(), PLY_ASCII, NULL, 0, NULL)) == NULL)
	{
		cout << "\nCannot open file for writing!!";
		exit(0);
	}
	
	// mark all darts defining infinite vertex
	int infiniteVertexMark = lcc.get_new_mark();
	
	if (infiniteVertexMark == -1)
		exit(0);

	// Marking all darts associated with infinite vertices
	markInfiniteVertexDart(dartToInfiniteVertex, lcc, infiniteVertexMark);
	
	// count number of vertices and faces in LCC
	size_t nVertices = 0, nFaces = 0;
	for (LCCWithIntInfo::One_dart_per_cell_range<0>::iterator pointCountIter = lcc.one_dart_per_cell<0>().begin(), pointCountIterEnd = lcc.one_dart_per_cell<0>().end(); pointCountIter != pointCountIterEnd; pointCountIter++)
		if (!isInfinite(pointCountIter, lcc, infiniteVertexMark, 0))
			nVertices++;

	for (LCCWithIntInfo::One_dart_per_cell_range<2>::iterator faceCountIter = lcc.one_dart_per_cell<2>().begin(), faceCountIterEnd = lcc.one_dart_per_cell<2>().end(); faceCountIter != faceCountIterEnd; faceCountIter++)
		if(!isInfinite(faceCountIter, lcc, infiniteVertexMark, 2))
			nFaces++;
	
	ply_add_element(lccOutputPLY, "vertex", nVertices);
	ply_add_scalar_property(lccOutputPLY, "x", PLY_FLOAT);
	ply_add_scalar_property(lccOutputPLY, "y", PLY_FLOAT);
	ply_add_scalar_property(lccOutputPLY, "z", PLY_FLOAT);

	ply_add_element(lccOutputPLY, "face", nFaces);
	ply_add_list_property(lccOutputPLY, "vertex_indices", PLY_UCHAR, PLY_INT32);

	if (!ply_write_header(lccOutputPLY))
	{
		cout << "Header cannot be written!!";
		exit(0);
	}

	// write vertices
	size_t pointId = 0;
	for (LCCWithIntInfo::One_dart_per_cell_range<0>::iterator pointIter = lcc.one_dart_per_cell<0>().begin(), pointIterEnd = lcc.one_dart_per_cell<0>().end(); pointIter != pointIterEnd; pointIter++)
	{
		if (!isInfinite(pointIter, lcc, infiniteVertexMark, 0))
		{
			CGALPoint pt = lcc.point(pointIter); 
			ply_write(lccOutputPLY, pt.x());
			ply_write(lccOutputPLY, pt.y());
			ply_write(lccOutputPLY, pt.z());
			lcc.info<0>(pointIter) = pointId++;
		}
	}
	
        // write polygons	
	size_t nFiniteFaces = 0, nInfiniteFaces = 0;
	for (LCCWithIntInfo::One_dart_per_cell_range<2>::iterator faceIter = lcc.one_dart_per_cell<2>().begin(), faceIterEnd = lcc.one_dart_per_cell<2>().end(); faceIter != faceIterEnd; faceIter++)
	{
		if (!isInfinite(faceIter, lcc, infiniteVertexMark, 2))
		{
			ply_write(lccOutputPLY, 3);
			for (LCCWithIntInfo::One_dart_per_incident_cell_range<0, 2>::iterator pointInFaceIter = lcc.one_dart_per_incident_cell<0, 2>(faceIter).begin(), pointInFaceIterEnd = lcc.one_dart_per_incident_cell<0, 2>(faceIter).end(); pointInFaceIter != pointInFaceIterEnd; pointInFaceIter++)
				ply_write(lccOutputPLY, lcc.info<0>(pointInFaceIter)); 
			nFiniteFaces++;
		}
		else
			nInfiniteFaces++;
	}

	lcc.free_mark(infiniteVertexMark);
	ply_close(lccOutputPLY);			
}



/*! \fn void CDTGenerator::computeDelaunayTetrahedralization()
    \brief Computes Delaunay tetrahedralization.
*/
void CDTGenerator::computeDelaunayTetrahedralization()
{	
	vector <pair <CGALPoint, DartHandle> > lccVertexVector;

	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
		lccVertexVector.push_back(pair<CGALPoint, DartHandle>(plc.point(pIter), pIter));

	DT.insert(lccVertexVector.begin(), lccVertexVector.end());

//	cout << "\nDelaunay tetrahedralization computed!!";

	LCCWithIntInfo DTLCC;
	string fileName("../../data/delaunay.ply");

	LCCWithIntInfo::Dart_handle dartToInfiniteVertex = import_from_triangulation_3(DTLCC, DT);
	
//	writePLYOutput(dartToInfiniteVertex, DTLCC, fileName);
}


/*! \fn void CDTGenerator::formMissingSegmentsQueue(vector<DartHandle> &missingSegmentQueue)
    \brief Collects missing constraint segments in a queue

    Determines which constraint segments are missing from Delaunay tetrahedralization of vertices of PLC.

    \param missingSegmentQueue [Out] Contains the _DartHandle_ to the missing constraint segments after execution.
*/
void CDTGenerator::formMissingSegmentsQueue(vector<DartHandle> &missingSegmentQueue)
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

	cout << "\nTotal number of missing constraint segments:" << missingSegmentQueue.size() << endl;


	return;
}


/*! \fn unsigned int CDTGenerator::computeCircumradius(CGALPoint &A, CGALPoint &B, CGALPoint &encroachingCandidate)
    \brief Computes circumradius of circle defined by input points.

    \param [in] A Endpoint 1 of constraint segment
    \param [in] B Endpoint 2 of constraint segment
    \param [in] encroachingCandidate Candidate for reference point
*/
unsigned int CDTGenerator::computeCircumradius(CGALPoint &A, CGALPoint &B, CGALPoint &encroachingCandidate)
{
	// computing circumradius of a triangle:
	
	float a, b, c;
	float circumradius;
	a = sqrt(pow((A.x() - encroachingCandidate.x()), 2)+ pow((A.y() - encroachingCandidate.y()), 2) + pow((A.z() - encroachingCandidate.z()), 2));
	b = sqrt(pow((B.x() - encroachingCandidate.x()), 2)+ pow((B.y() - encroachingCandidate.y()), 2) + pow((B.z() - encroachingCandidate.z()), 2));
	c = sqrt(pow((A.x() - B.x()), 2)+ pow((A.y() - B.y()), 2) + pow((A.z() - B.z()), 2));

	return circumradius = (a * b * c) / sqrt((a + b + c) * (b + c - a) * (c + a - b) * (a + b - c));
}


/*! \fn void CDTGenerator::computeReferencePoint(CGALPoint *refPoint, DartHandle missingSegmentHandle)
    \brief Computes the reference point for segment splitting

    Computes reference point which will be used for determing the location of _steiner point_ for the input missing constraint segment.

    \param [out] refPoint Pointer to the computed reference point.
    \param [in] missingSegmentHandle Dart handle of missing constraint segment.
*/
void CDTGenerator::computeReferencePoint(CGALPoint *refPoint, DartHandle missingSegmentHandle)
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


/*! \fn float CDTGenerator::dotProduct(DartHandle segment1Handle, DartHandle segment2Handle)
    \brief Computes dot product of vectors represented by input segments.

    \param [in] segment1Handle DartHandle for the first segment
    \param [in] segment2Handle DartHandle for the second segment
*/
float CDTGenerator::dotProduct(DartHandle segment1Handle, DartHandle segment2Handle)
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


/*! \fn float CDTGenerator::vectorMagnitude(DartHandle inputSegmentHandle)
    \brief Computes magnitue of input vector

    \param [in] inputSegmentHandle DartHandle of input segment
*/
float CDTGenerator::vectorMagnitude(DartHandle inputSegmentHandle)
{
	CGALPoint segmentVertices[2];

	segmentVertices[0] = plc.point(inputSegmentHandle);
	segmentVertices[1] = plc.point(plc.beta(inputSegmentHandle, 1));	

	CGALPoint vector = CGALPoint(segmentVertices[0].x() - segmentVertices[1].x(), segmentVertices[0].y() - segmentVertices[1].y(), segmentVertices[0].z() - segmentVertices[1].z());

	float vectorMagnitude = sqrtf(powf(vector.x(), 2.0) + powf(vector.y(), 2.0) + powf(vector.z(), 2.0));

	return vectorMagnitude;
}


/*! \fn float CDTGenerator::computeAngleBetweenSegments(DartHandle segment1Handle, DartHandle segment2Handle)
    \brief Computes angle between vectors represented by input segments.

    \param [in] segment1Handle DartHandle for first input segment
    \param [in] segment2Handle DartHandle for second input segment
*/
float CDTGenerator::computeAngleBetweenSegments(DartHandle segment1Handle, DartHandle segment2Handle)
{

	float angle = acosf(dotProduct(segment1Handle, segment2Handle) / (vectorMagnitude(segment1Handle) * vectorMagnitude(segment2Handle)));

	return (angle * 180.0f / Pi); // conversion to degrees
}


/*! \fn bool CDTGenerator::isVertexAcute(DartHandle inputPointHandle)
    \brief Tests whether input vertex is acute

    A vertex is _not_ acute if there exists _atleast_ one pair of _incident_ segments which have angle _greater_ than 90 degrees.    
    \param [in] inputPointHandle DartHandle of input vertex
*/
bool CDTGenerator::isVertexAcute(DartHandle inputPointHandle)
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


/*! \fn unsigned int CDTGenerator::determineSegmentType(DartHandle missingSegmentHandle)
    \brief Determines the input segment type

    A segment is of _type-1_ if _none_ of its endpoints are _acute_. If _exactly_ one endpoint is acute than segment is of _type-2_.

    \param [in] DartHandle for input constraint segment
*/
unsigned int CDTGenerator::determineSegmentType(DartHandle missingSegmentHandle)
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


/*! \fn float CDTGenerator::computeSegmentLength(CGALPoint &A, CGALPoint &B)
    \brief Computes Euler length of the segment represented by _A_ and _B_.

    \param [in] A First endpoint of the segment
    \param [in] B Second endpoint of the segment
*/
float CDTGenerator::computeSegmentLength(CGALPoint &A, CGALPoint &B)
{
	float sLength;
	sLength = sqrt(pow(A.x() - B.x(), 2) + pow(A.y() - B.y(), 2) + pow(A.z() - B.z(), 2));
	return sLength;
}


/*! \fn void CDTGenerator::updatePLCAndDT(CGALPoint &v, DartHandle missingSegmentHandle)
   /brief Updates PLC and Delaunay triangulation after insertion of a new vertex into a constraint segment.

   /param [in] v New vertex to be inserted into a constraint segment of PLC
   /param [ih] missignSegmentHandle DartHandle of constraint segment which will be split after inserting new vertex  
*/
void CDTGenerator::updatePLCAndDT(CGALPoint &v, DartHandle missingSegmentHandle)
{

	// update PLC
	plc.insert_point_in_cell<1>(missingSegmentHandle, v);
	// update DT
	computeDelaunayTetrahedralization(); 
}


/*! \fn void CDTGenerator::splitMissingSegment(DartHandle missingSegmentHandle)
    \brief Splits the missing constraint segment.

    A constraint segment is split(by insertion of a new vertex) by applying segment splitting rules. Purpose of splitting is to make resulting subsegments _strongly Delaunay_.

    \param [in] missingSegmentHandle DartHandle of missing constraint segment
*/
void CDTGenerator::splitMissingSegment(DartHandle missingSegmentHandle)
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


/*! \fn void CDTGenerator::recoverConstraintSegments()
    \brief Top-level routine for recovering constraint segments.

    Missing constriant segments of PLC are recovered by spliting the segments in such a way so that resulting sub-segments become _strongly Delaunay_.
*/
void CDTGenerator::recoverConstraintSegments()
{
	vector<DartHandle> missingSegmentQueue;
	DartHandle missingSegment;

	formMissingSegmentsQueue(missingSegmentQueue);
	int i = 0;
	while (missingSegmentQueue.size() != 0)
	{
		cout << "Segment recovery iteration: #" << i << endl;
		missingSegment = missingSegmentQueue.back();
		missingSegmentQueue.pop_back();
		splitMissingSegment(missingSegment);
		formMissingSegmentsQueue(missingSegmentQueue); // update missingSegmentQueue
		i++;
	}

	return;
}


/*! \fn void CDTGenerator::generate()
    \brief Public interface for CDTGenerator class.
 */
void CDTGenerator::generate()
{
	readPLCInput();
	computeDelaunayTetrahedralization();
	recoverConstraintSegments();
}


