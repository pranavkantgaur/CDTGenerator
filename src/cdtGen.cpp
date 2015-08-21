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
bool CDTGenerator::areGeometricallySameSegments(DartHandle d1, DartHandle d2, LCC& lcc)
{
	if (lcc.point(d1) == lcc.point(lcc.beta(d2, 1)))
		if (lcc.point(lcc.beta(d1, 1)) == lcc.point(d2))
			return true;
	if (lcc.point(d1) == lcc.point(d2))
		if (lcc.point(lcc.beta(d1, 1)) == lcc.point(lcc.beta(d2, 1)))
			return true;
	else
		return false;
} 


/*! \fn bool CDTGenerator::areGeometricallySameSegmentsWithDartInfo(DartHandle d1, DartHandle d2)
    \brief Tests whether segments represented by d1 and d2 are _geometrically_ same.
    \param [in] d1 DartHandle for first segment
    \param [in] d2 DartHandle for second segment
*/
bool CDTGenerator::areGeometricallySameSegmentsWithDartInfo(LCCWithDartInfo::Dart_handle d1, LCCWithDartInfo::Dart_handle d2, LCCWithDartInfo &lcc)
{
	if (lcc.point(d1) == lcc.point(lcc.beta(d2, 1)))
		if (lcc.point(lcc.beta(d1, 1)) == lcc.point(d2))
			return true;
	return false;
} 


/*! \fn void CDTGenerator::sew2CellsFromEdge(LCC &lcc)
 *  \brief Sews 2-cells sharing common edge.
 *  \param [in, out] lcc Linear cell complex containing 2-cells. 
 */
void CDTGenerator::sew2CellsFromEdge(LCC &lcc)
{

	vector<pair<DartHandle, DartHandle> > twoCellsToBeSewed;
	int sewedMark = lcc.get_new_mark();

	if (sewedMark == -1)
	{
		cout << "\nNo free mark available!!";
		exit(0);
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
					if (areGeometricallySameSegments(segIter1, segIter2, lcc)) // checks the geometry of segments
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
	size_t k = 0;
	for (vector<pair<DartHandle, DartHandle> >::iterator dIter = twoCellsToBeSewed.begin(), dIterEnd = twoCellsToBeSewed.end(); dIter != dIterEnd; dIter++)
		if (lcc.is_sewable<2>(dIter->first, dIter->second))
		{
			lcc.sew<2>(dIter->first, dIter->second);
			k++;	
		}

	lcc.free_mark(sewedMark);
}


/*! \fn void CDTGenerator::sew2CellsWithDartInfoFromEdge(LCCWithDartInfo &lcc)
 *  \brief Sews 2-cells of LCC with dart info sharing common edge.
 *  \param [in, out] Linear cell complex containing 2-cells.  
 */
void CDTGenerator::sew2CellsWithDartInfoFromEdge(LCCWithDartInfo &lcc)
{

	vector<pair<LCCWithDartInfo::Dart_handle, LCCWithDartInfo::Dart_handle> > twoCellsToBeSewed;
	int sewedMark = lcc.get_new_mark();

	if (sewedMark == -1)
	{
		cout << "\nNo free mark available!!";
		exit(0);
	}	

	// sew facets sharing edge
	for (LCCWithDartInfo::Dart_range::iterator segIter1 = lcc.darts().begin(), segIterEnd1 = lcc.darts().end(); segIter1 != segIterEnd1; segIter1++)
	{
		if (!lcc.is_marked(segIter1, sewedMark)) // not sewed till now
		{
			for (LCCWithDartInfo::Dart_range::iterator segIter2 = lcc.darts().begin(), segIterEnd2 = lcc.darts().end(); segIter2 != segIterEnd2; segIter2++)
			{
				if (!lcc.is_marked(segIter2, sewedMark) && lcc.is_sewable<2>(segIter1, segIter2))
				{
					if (areGeometricallySameSegmentsWithDartInfo(segIter1, segIter2, lcc)) // checks the geometry of segments
					{
						lcc.mark(segIter1, sewedMark);
						lcc.mark(segIter2, sewedMark);
						twoCellsToBeSewed.push_back(pair<LCCWithDartInfo::Dart_handle, LCCWithDartInfo::Dart_handle>(segIter1, segIter2));
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
	size_t k = 0;
	for (vector<pair<LCCWithDartInfo::Dart_handle, LCCWithDartInfo::Dart_handle> >::iterator dIter = twoCellsToBeSewed.begin(), dIterEnd = twoCellsToBeSewed.end(); dIter != dIterEnd; dIter++)
		if (lcc.is_sewable<2>(dIter->first, dIter->second))
		{
			lcc.sew<2>(dIter->first, dIter->second);
			k++;	
		}
	lcc.free_mark(sewedMark);
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
	size_t vertexIds[3];
 
	CGALPoint trianglePoints[3];

	for (size_t n = 0, m = plcFaceVector.size(); n < m; n++)
	{
		for (size_t k = 0; k < 3; k++)
			vertexIds[k] = plcFaceVector[n].pointIds[k];
	
		for (size_t i = 0; i < 3; i++)
			trianglePoints[i] = CGALPoint(plcVertexVector[vertexIds[i]].x(), plcVertexVector[vertexIds[i]].y(), plcVertexVector[vertexIds[i]].z());
			
		plc.make_triangle(trianglePoints[0], trianglePoints[1], trianglePoints[2]);
	}
	
	// sew facets sharing edge
	sew2CellsFromEdge(plc);
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


/*! void CDTGenerator::markInfinteVertexDart(LCC::Dart_handle d, LCC &lcc, int infiniteVertexMark)
*   \brief Marks all darts defining the infinite vertex in LCC.
*   \param [in] d Handle to the dart representing infinite vertex.
*   \param [in] lcc Linear cell complex to be marked.
*   \param [in] infiniteVertexMark Mark representing a dart defining the infinite vertex.
*/
void CDTGenerator::markInfiniteVertexDart(LCC::Dart_handle d, LCC &lcc, int infiniteVertexMark)
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


/*! \fn bool CDTGenerator::isInfinite(LCCWithIntInfo::Dart_handle adart, LCCWithIntInfo& lcc, int infiniteVertexMark, size_t cell_dimension)
 *  \brief Tests whether the given i-cell is infinite.
 *  \param [in] adart a dart handle to the i-cell.
 *  \param [in] lcc Linear cell complex to be checked.
 *  \param [in] infiniteVertexMark mark representing the infinite vertex.
 *  \param [in] cell_dimension dimension of i-cell.
 *  \return True if input i-cell is infinite.
 */
bool CDTGenerator::isInfinite(LCCWithIntInfo::Dart_handle adart, const LCCWithIntInfo& lcc, int infiniteVertexMark, size_t cell_dimension)
{
	if (infiniteVertexMark == INVALID_VALUE)
		return false;

	bool isInfinite = false;
	
	if (cell_dimension == 0)
	{
		if (lcc.is_marked(adart, infiniteVertexMark)) 
			return true;
	}

	if (cell_dimension == 2)
	{
		for (LCCWithIntInfo::Dart_of_orbit_const_range<1>::const_iterator pIter = lcc.darts_of_orbit<1>(adart).begin(), pIterEnd = lcc.darts_of_orbit<1>(adart).end(); pIter != pIterEnd; pIter++)
		{
			if (lcc.is_marked(pIter, infiniteVertexMark)) 
			{
				isInfinite = true;
				break;
			}
		}
	
	}

	if (cell_dimension == 3)
	{
		for (LCCWithIntInfo::One_dart_per_incident_cell_const_range<0, 3>::const_iterator pIter = lcc.one_dart_per_incident_cell<0, 3>(adart).begin(), pIterEnd = lcc.one_dart_per_incident_cell<0, 3>(adart).end(); pIter != pIterEnd; pIter++)
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



/*! \fn bool CDTGenerator::isInfinite(LCC::Dart_handle adart, LCC& lcc, int infiniteVertexMark, size_t cell_dimension)
 *  \brief Tests whether the given i-cell is infinite.
 *  \param [in] adart a dart handle to the i-cell.
 *  \param [in] lcc Linear cell complex to be checked.
 *  \param [in] infiniteVertexMark mark representing the infinite vertex.
 *  \param [in] cell_dimension dimension of i-cell.
 *  \return True if input i-cell is infinite.
*/
bool CDTGenerator::isInfinite(LCC::Dart_handle adart, const LCC& lcc, int infiniteVertexMark, size_t cell_dimension)
{

	if (infiniteVertexMark == INVALID_VALUE)
		return false;

	bool isInfinite = false;

	if (cell_dimension == 0)
	{
	
		if (lcc.is_marked(adart, infiniteVertexMark)) 
			return true;
	}

	if (cell_dimension == 2)
	{
		for (LCC::Dart_of_orbit_const_range<1>::const_iterator pIter = lcc.darts_of_orbit<1>(adart).begin(), pIterEnd = lcc.darts_of_orbit<1>(adart).end(); pIter != pIterEnd; pIter++)
		{
			if (lcc.is_marked(pIter, infiniteVertexMark)) 
			{
				isInfinite = true;
				break;
			}
		}
	
	}

	if (cell_dimension == 3)
	{
		for (LCC::One_dart_per_incident_cell_const_range<0, 3>::const_iterator pIter = lcc.one_dart_per_incident_cell<0, 3>(adart).begin(), pIterEnd = lcc.one_dart_per_incident_cell<0, 3>(adart).end(); pIter != pIterEnd; pIter++)
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


/*! \fn void CDTGenerator::copyLCCToLCCWithIntInfo()
 */
void CDTGenerator::copyLCCToLCCWithIntInfo(LCC &lcc, LCCWithIntInfo &lccWithIntInfo)
{

	CGALPoint p[4];
	size_t i;
	//copy tets
	for (LCC::One_dart_per_cell_range<3>::iterator cellHandle = lcc.one_dart_per_cell<3>().begin(), cellHandleEnd = lcc.one_dart_per_cell<3>().end(); cellHandle != cellHandleEnd; cellHandle++)
	{
		i = 0;
		for (LCC::One_dart_per_incident_cell_range<0, 3>::iterator pointHandle = lcc.one_dart_per_incident_cell<0, 3>(cellHandle).begin(), pointHandleEnd = lcc.one_dart_per_incident_cell<0, 3>(cellHandle).end(); pointHandle != pointHandleEnd; pointHandle++)
			p[i++] = lcc.point(pointHandle);
		
		lccWithIntInfo.make_tetrahedron(p[0], p[1], p[2], p[3]);
	}

}


/*! \fn void CDTGenerator::writePLYOutput(LCC::Dart_handle dartToInfiniteVertex, LCC &lcc, string fileName)
*   \brief Writes mesh represented as LCC to PLY file
*   \param [in] dartToInfiniteVertex Dart handle to the infinite vertex.
*   \param [in] lcc Linear cell complex.
*   \param [in] fileName Name of output PLY file.
*/
void CDTGenerator::writePLYOutput(LCC::Dart_handle dartToInfiniteVertex, LCC &lcc, string fileName)
{
	p_ply lccOutputPLY;
        
	if ((lccOutputPLY = ply_create(fileName.c_str(), PLY_ASCII, NULL, 0, NULL)) == NULL)
	{
		cout << "\nCannot open file for writing!!";
		exit(0);
	}

	LCCWithIntInfo lccWithIntInfo;
	copyLCCToLCCWithIntInfo(lcc, lccWithIntInfo);


	int infiniteVertexMark;

	if (dartToInfiniteVertex == NULL)
		infiniteVertexMark = INVALID_VALUE;
	else
	{
		// mark all darts defining infinite vertex
		infiniteVertexMark = lccWithIntInfo.get_new_mark();
	}
	
	if (infiniteVertexMark == -1)
		exit(0);
	
	// Marking all darts associated with infinite vertices
	markInfiniteVertexDart(NULL, lccWithIntInfo, infiniteVertexMark);
	
	// count number of vertices and faces in LCC
	size_t nVertices = 0, nTets = 0;
	for (LCCWithIntInfo::One_dart_per_cell_range<0>::iterator pointCountIter = lccWithIntInfo.one_dart_per_cell<0>().begin(), pointCountIterEnd = lccWithIntInfo.one_dart_per_cell<0>().end(); pointCountIter != pointCountIterEnd; pointCountIter++)
		if (!isInfinite(pointCountIter, lccWithIntInfo, infiniteVertexMark, 0))
			nVertices++;

	for (LCCWithIntInfo::One_dart_per_cell_range<3>::iterator faceCountIter = lccWithIntInfo.one_dart_per_cell<3>().begin(), faceCountIterEnd = lccWithIntInfo.one_dart_per_cell<3>().end(); faceCountIter != faceCountIterEnd; faceCountIter++)
		if(!isInfinite(faceCountIter, lccWithIntInfo, infiniteVertexMark, 3))
			nTets++;
	
	ply_add_element(lccOutputPLY, "vertex", nVertices);
	ply_add_scalar_property(lccOutputPLY, "x", PLY_FLOAT);
	ply_add_scalar_property(lccOutputPLY, "y", PLY_FLOAT);
	ply_add_scalar_property(lccOutputPLY, "z", PLY_FLOAT);

	ply_add_element(lccOutputPLY, "face", nTets*4);
	ply_add_list_property(lccOutputPLY, "vertex_indices", PLY_UCHAR, PLY_INT32);

	if (!ply_write_header(lccOutputPLY))
	{
		cout << "Header cannot be written!!";
		exit(0);
	}
	
//	cout << "#### Writing vertices..." << endl;
	// write vertices
	size_t pointId = 0;
	for (LCCWithIntInfo::One_dart_per_cell_range<0>::iterator pointIter = lccWithIntInfo.one_dart_per_cell<0>().begin(), pointIterEnd = lccWithIntInfo.one_dart_per_cell<0>().end(); pointIter != pointIterEnd; pointIter++)
	{
		if (!isInfinite(pointIter, lcc, infiniteVertexMark, 0))
		{
			CGALPoint pt = lccWithIntInfo.point(pointIter); 
			ply_write(lccOutputPLY, pt.x());
			ply_write(lccOutputPLY, pt.y());
			ply_write(lccOutputPLY, pt.z());
			lccWithIntInfo.info<0>(pointIter) = pointId++;
		}
	}
	//sew tets
	lccWithIntInfo.sew3_same_facets();

	// write tetrahedrons	
//	cout << "#### Writing polygons..." << endl;
	LCCWithIntInfo::Dart_handle pts[4];
	size_t i;
 
	size_t nFiniteTets = 0, nInfiniteTets = 0;
	for (LCCWithIntInfo::One_dart_per_cell_range<3>::iterator tetIter = lccWithIntInfo.one_dart_per_cell<3>().begin(), tetIterEnd = lccWithIntInfo.one_dart_per_cell<3>().end(); tetIter != tetIterEnd; tetIter++)
	{
		if (!isInfinite(tetIter, lccWithIntInfo, infiniteVertexMark, 3))
		{
			i = 0;
 
			for (LCCWithIntInfo::One_dart_per_incident_cell_range<0, 3>::iterator pIter = lccWithIntInfo.one_dart_per_incident_cell<0, 3>(tetIter).begin(), pIterEnd = lccWithIntInfo.one_dart_per_incident_cell<0, 3>(tetIter).end(); pIter != pIterEnd; pIter++)
				pts[i++] = pIter;
			// t1	
			ply_write(lccOutputPLY, 3);	
			ply_write(lccOutputPLY, lccWithIntInfo.info<0>(pts[0]));
			ply_write(lccOutputPLY, lccWithIntInfo.info<0>(pts[1]));
			ply_write(lccOutputPLY, lccWithIntInfo.info<0>(pts[2]));
		
			// t2
			ply_write(lccOutputPLY, 3);
			ply_write(lccOutputPLY, lccWithIntInfo.info<0>(pts[1]));
			ply_write(lccOutputPLY, lccWithIntInfo.info<0>(pts[0]));
			ply_write(lccOutputPLY, lccWithIntInfo.info<0>(pts[3]));
		
			// t3 	
			ply_write(lccOutputPLY, 3);	
			ply_write(lccOutputPLY, lccWithIntInfo.info<0>(pts[1]));
			ply_write(lccOutputPLY, lccWithIntInfo.info<0>(pts[3]));
			ply_write(lccOutputPLY, lccWithIntInfo.info<0>(pts[2]));

			// t4
			ply_write(lccOutputPLY, 3);	
			ply_write(lccOutputPLY, lccWithIntInfo.info<0>(pts[3]));
			ply_write(lccOutputPLY, lccWithIntInfo.info<0>(pts[0]));
			ply_write(lccOutputPLY, lccWithIntInfo.info<0>(pts[2]));
		}
	}

	
//	cout << "Output file written successfully!" << endl;
	lccWithIntInfo.free_mark(infiniteVertexMark);
	ply_close(lccOutputPLY);			

}


/*! \fn void CDTGenerator::writePLYOutput(LCCWithIntInfo::Dart_handle dartToInfiniteVertex, LCCWithIntInfo &lcc, string fileName)
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
	size_t nVertices = 0, nTets = 0;
	for (LCCWithIntInfo::One_dart_per_cell_range<0>::iterator pointCountIter = lcc.one_dart_per_cell<0>().begin(), pointCountIterEnd = lcc.one_dart_per_cell<0>().end(); pointCountIter != pointCountIterEnd; pointCountIter++)
		if (!isInfinite(pointCountIter, lcc, infiniteVertexMark, 0))
			nVertices++;

	for (LCCWithIntInfo::One_dart_per_cell_range<3>::iterator faceCountIter = lcc.one_dart_per_cell<3>().begin(), faceCountIterEnd = lcc.one_dart_per_cell<3>().end(); faceCountIter != faceCountIterEnd; faceCountIter++)
		if(!isInfinite(faceCountIter, lcc, infiniteVertexMark, 3))
			nTets++;
	
	ply_add_element(lccOutputPLY, "vertex", nVertices);
	ply_add_scalar_property(lccOutputPLY, "x", PLY_FLOAT);
	ply_add_scalar_property(lccOutputPLY, "y", PLY_FLOAT);
	ply_add_scalar_property(lccOutputPLY, "z", PLY_FLOAT);

	ply_add_element(lccOutputPLY, "face", nTets*4);
	ply_add_list_property(lccOutputPLY, "vertex_indices", PLY_UCHAR, PLY_INT32);

	if (!ply_write_header(lccOutputPLY))
	{
		cout << "Header cannot be written!!";
		exit(0);
	}
	
//	cout << "#### Writing vertices..." << endl;
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
	
	// write tetrahedrons	
//	cout << "#### Writing polygons..." << endl;
	LCCWithIntInfo::Dart_handle pts[4];
	size_t i;
 
	size_t nFiniteTets = 0, nInfiniteTets = 0;
	for (LCCWithIntInfo::One_dart_per_cell_range<3>::iterator tetIter = lcc.one_dart_per_cell<3>().begin(), tetIterEnd = lcc.one_dart_per_cell<3>().end(); tetIter != tetIterEnd; tetIter++)
	{
		if (!isInfinite(tetIter, lcc, infiniteVertexMark, 3))
		{
			i = 0;
 
			for (LCCWithIntInfo::One_dart_per_incident_cell_range<0, 3>::iterator pIter = lcc.one_dart_per_incident_cell<0, 3>(tetIter).begin(), pIterEnd = lcc.one_dart_per_incident_cell<0, 3>(tetIter).end(); pIter != pIterEnd; pIter++)
				pts[i++] = pIter;
			// t1	
			ply_write(lccOutputPLY, 3);	
			ply_write(lccOutputPLY, lcc.info<0>(pts[0]));
			ply_write(lccOutputPLY, lcc.info<0>(pts[1]));
			ply_write(lccOutputPLY, lcc.info<0>(pts[2]));
		
			// t2
			ply_write(lccOutputPLY, 3);
			ply_write(lccOutputPLY, lcc.info<0>(pts[1]));
			ply_write(lccOutputPLY, lcc.info<0>(pts[0]));
			ply_write(lccOutputPLY, lcc.info<0>(pts[3]));
		
			// t3 	
			ply_write(lccOutputPLY, 3);	
			ply_write(lccOutputPLY, lcc.info<0>(pts[1]));
			ply_write(lccOutputPLY, lcc.info<0>(pts[3]));
			ply_write(lccOutputPLY, lcc.info<0>(pts[2]));

			// t4
			ply_write(lccOutputPLY, 3);	
			ply_write(lccOutputPLY, lcc.info<0>(pts[3]));
			ply_write(lccOutputPLY, lcc.info<0>(pts[0]));
			ply_write(lccOutputPLY, lcc.info<0>(pts[2]));
		}
	}

	
//	cout << "Output file written successfully!" << endl;
	lcc.free_mark(infiniteVertexMark);
	ply_close(lccOutputPLY);			

}



/*! \fn void CDTGenerator::computeDelaunayTetrahedralization()
    \brief Computes Delaunay tetrahedralization.
*/
void CDTGenerator::computeDelaunayTetrahedralization(int missingSegmentQueueSize)
{	
	vector <pair <CGALPoint, DartHandle> > lccVertexVector;

	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
		lccVertexVector.push_back(pair<CGALPoint, DartHandle>(plc.point(pIter), pIter));

	DT.insert(lccVertexVector.begin(), lccVertexVector.end());

//	cout << "\nDelaunay tetrahedralization computed!!";

	LCCWithIntInfo DTLCC;
	string fileName("../../data/delaunay.ply");

	LCCWithIntInfo::Dart_handle dartToInfiniteVertex = import_from_triangulation_3(DTLCC, DT);
	
//	cout << "Writing Delaunay triangulation to output file..." << endl;	
	if (missingSegmentQueueSize == 0)
		writePLYOutput(dartToInfiniteVertex, DTLCC, fileName);
}


/*! \fn void CDTGenerator::formMissingSegmentsQueue(vector<DartHandle> &missingSegmentQueue)
    \brief Collects missing constraint segments in a queue

    Determines which constraint segments are missing from Delaunay tetrahedralization of vertices of PLC.

    \param missingSegmentQueue [Out] Contains the _DartHandle_ to the missing constraint segments after execution.
*/
void CDTGenerator::formMissingSegmentsQueue()
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


/*! \fn size_t CDTGenerator::computeCircumradius(CGALPoint &A, CGALPoint &B, CGALPoint &encroachingCandidate)
    \brief Computes circumradius of circle defined by input points.

    \param [in] A Endpoint 1 of constraint segment
    \param [in] B Endpoint 2 of constraint segment
    \param [in] encroachingCandidate Candidate for reference point
*/
size_t CDTGenerator::computeCircumradius(CGALPoint &A, CGALPoint &B, CGALPoint &encroachingCandidate)
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

	CGALPoint vector1 = CGALPoint(fabs(segment1Vertex[0].x() - segment1Vertex[1].x()), fabs(segment1Vertex[0].y() - segment1Vertex[1].y()), fabs(segment1Vertex[0].z() - segment1Vertex[1].z()));
	CGALPoint vector2 = CGALPoint(fabs(segment2Vertex[0].x() - segment2Vertex[1].x()), fabs(segment2Vertex[0].y() - segment2Vertex[1].y()), fabs(segment2Vertex[0].z() - segment2Vertex[1].z()));

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

	// determine all segments incident on input vertex 
	for (LCC::One_dart_per_incident_cell_range<1, 0>::iterator incidentSegmentIter = plc.one_dart_per_incident_cell<1, 0>(inputPointHandle).begin(), incidentSegmentIterEnd = plc.one_dart_per_incident_cell<1, 0>(inputPointHandle).end(); incidentSegmentIter != incidentSegmentIterEnd;incidentSegmentIter++)
		incidentOnInputPoint.push_back(incidentSegmentIter);

	// Compute angle between all possible pairs(TODO: NAIVE SOLUTION)
	for (vector<DartHandle>::iterator segIter1 = incidentOnInputPoint.begin(); segIter1 != incidentOnInputPoint.end(); segIter1++) 
		for (vector<DartHandle>::iterator segIter2 = incidentOnInputPoint.begin(); *segIter1 != *segIter2 && segIter2 != incidentOnInputPoint.end(); segIter2++)
			if (computeAngleBetweenSegments(*segIter1, *segIter2) < 90.0f)
				return true;

	return false; 
}


/*! \fn size_t CDTGenerator::determineSegmentType(DartHandle missingSegmentHandle)
    \brief Determines the input segment type

    A segment is of _type-1_ if _none_ of its endpoints are _acute_. If _exactly_ one endpoint is acute than segment is of _type-2_.

    \param [in] DartHandle for input constraint segment
*/
size_t CDTGenerator::determineSegmentType(DartHandle missingSegmentHandle)
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

	// TODO: update PLC 
	//// remove all 2-cells containing this edge
	//// for each removed triangle (A, B, C): add triangles as: (A, B, C) + v => (A, B, v) + (A, v, C)
	//// stitch triangles sharing common edge
	
/*	cout << "DEBUG(Before plc update) !!" << endl;
	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
		cout << plc.point(pIter) << endl;
 	cout << "ENDS!!" << endl;
*/	
	vector<LCC::Dart_handle> incidentFacets;
	for (LCC::One_dart_per_incident_cell_range<2, 1>::iterator fIter = plc.one_dart_per_incident_cell<2, 1>(missingSegmentHandle).begin(), fIterEnd = plc.one_dart_per_incident_cell<2, 1>(missingSegmentHandle).end(); fIter != fIterEnd; fIter++)
		incidentFacets.push_back(fIter);

	CGALPoint A, B, C;
	for (vector<LCC::Dart_handle>::iterator facetIter = incidentFacets.begin(), facetIterEnd = incidentFacets.end(); facetIter != facetIterEnd; facetIter++)
	{
		for (LCC::One_dart_per_incident_cell_range<1, 2>::iterator sIter = plc.one_dart_per_incident_cell<1, 2>(*facetIter).begin(), sIterEnd = plc.one_dart_per_incident_cell<1, 2>(*facetIter).end(); sIter != sIterEnd; sIter++)	
		{

			if (areGeometricallySameSegments(sIter, missingSegmentHandle, plc))
			{
				A = plc.point(plc.beta<1, 1>(sIter));
				B = plc.point(sIter);
				C = plc.point(plc.beta<1>(sIter));
				plc.make_triangle(A, B, v);
				plc.make_triangle(A, v, C);
			}
			
		}
			remove_cell<LCC, 2>(plc, *facetIter); 
			sew2CellsFromEdge(plc);

	/*		cout << "Geometrically identical segment found!!" << endl;
			cout << "v: " << endl;
			cout << "A: " << A << endl;
			cout << "B: " << B << endl;
			cout << "C: " << C << endl;*/
	}	
	
	/*cout << "DEBUG(After plc update) !!" << endl;
	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
		cout << plc.point(pIter) << endl;
 	cout << "ENDS!!" << endl;
	*/
	
	// update DT
	computeDelaunayTetrahedralization(missingSegmentQueue.size()); 
}


/*! \fn void CDTGenerator::splitMissingSegment(DartHandle missingSegmentHandle)
    \brief Splits the missing constraint segment.

    A constraint segment is split(by insertion of a new vertex) by applying segment splitting rules. Purpose of splitting is to make resulting subsegments _strongly Delaunay_.

    \param [in] missingSegmentHandle DartHandle of missing constraint segment
*/
void CDTGenerator::splitMissingSegment(DartHandle missingSegmentHandle)
{
	cout << "Inside segment splitting function!!" << endl;
	CGALPoint vb, refPoint;
	CGALPoint sphereCenter;
	float sphereRadius;
	size_t segmentType;

	// determine segment type
	segmentType = determineSegmentType(missingSegmentHandle);

	// compute reference point
	computeReferencePoint(&refPoint, missingSegmentHandle);
	
	// endpoints of missing segment	
	CGALPoint A = plc.point(missingSegmentHandle);
	CGALPoint B = plc.point(plc.beta(missingSegmentHandle, 1));
	
//	cout << "Point A: " << A << endl;
//	cout << "Point B: " << B << endl;

	float AP, PB, AB;

	CGALPoint v; // steiner point

	if (segmentType == 1)
	{
		AP = computeSegmentLength(A, refPoint);
		AB = computeSegmentLength(A, B);
		PB = computeSegmentLength(refPoint, B);
		
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

		// Compute coordinates of steiner point:
		CGALSphericalPoint sphericalSphereCenter = CGALSphericalPoint(sphereCenter.x(), sphereCenter.y(), sphereCenter.z());
		CGALSphericalSphere s = CGALSphericalSphere(sphericalSphereCenter, pow(sphereRadius, 2));

		CGALSphericalPoint p1 = CGALSphericalPoint(A.x(), A.y(), A.z());
		CGALSphericalPoint p2 = CGALSphericalPoint(B.x(), B.y(), B.z());

		CGALSphericalSegment seg(p1, p2);
		CGALSphericalLineArc lineArc(seg);
		vector<Object> intersections; 
//		cout << "Case 1!!" << endl;		
		intersection(lineArc, s, back_inserter(intersections));
//		cout << "Case 1 ends!!" << endl;
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
		// identify the acute vertex out of A,B or one of their parents
		DartHandle acuteParentHandle; 
		DartHandle ApB[2]; 
		float vbLength;
		DartHandle AHandle = missingSegmentHandle;
		DartHandle BHandle = plc.beta(missingSegmentHandle, 1);

		if (isVertexAcute(AHandle) && !isVertexAcute(BHandle)) // if A is acute but B is not
		{
			acuteParentHandle = AHandle;
//			cout << "A is acute!!" << endl;
		}
		else if (isVertexAcute(BHandle) && !isVertexAcute(AHandle)) // if B is acute but A is not
		{
//			cout << "B is acute!!" << endl;
			acuteParentHandle = BHandle;
			BHandle = AHandle;
			AHandle = acuteParentHandle;
		}
		
//		cout << "Point A: " << plc.point(AHandle) << endl;
//	        cout << "Point B: " << plc.point(BHandle) << endl;	
	
		// Segment calculations		
		ApB[0] = acuteParentHandle;
		ApB[1] = BHandle;
			
		CGALSphericalPoint p1(plc.point(ApB[0]).x(), plc.point(ApB[0]).y(), plc.point(ApB[0]).z());
		CGALSphericalPoint p2(plc.point(ApB[1]).x(), plc.point(ApB[1]).y(), plc.point(ApB[1]).z());
		
		CGALSphericalSegment seg(p1, p2);
		CGALSphericalLineArc lineArc(seg);
		
		/// Sphere calculations
		CGALSphericalPoint acuteParent(plc.point(acuteParentHandle).x(), plc.point(acuteParentHandle).y(), plc.point(acuteParentHandle).z());
		
		CGALPoint acuteParentLinearField(plc.point(acuteParentHandle).x(), plc.point(acuteParentHandle).y(), plc.point(acuteParentHandle).z());
	
		float ARefPointLength = computeSegmentLength(acuteParentLinearField, refPoint);

		CGALSphericalSphere s(acuteParent, pow(ARefPointLength, 2));

		vector<Object> intersections;
//		cout << "Case 2!!" << endl;
//		cout << "Endpoint 1: " << p1 << endl;
//		cout << "Endpoint 2: " << p2 << endl;
		
		intersection(lineArc, s, back_inserter(intersections));
//		cout << "Case 2 ends!!" << endl;
	
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
		
		vbLength = computeSegmentLength(v, plc.point(BHandle)); 

		if (vbLength < vrefpointLength) // v was rejected
		{
			CGALSphericalPoint sphereCenter(acuteParent);
			size_t avLength = computeSegmentLength(plc.point(AHandle), v);
			if (vrefpointLength < 0.5 * avLength)
			{
					CGALPoint temp(acuteParentLinearField.x(), acuteParentLinearField.y(), acuteParentLinearField.z());
					acuteparentALength = computeSegmentLength(temp, plc.point(AHandle));
					sphereRadius = acuteparentALength + avLength - vrefpointLength;	
			}
			else
				sphereRadius = acuteparentALength + 0.5 * avLength;
		        
			s = CGALSphericalSphere(sphereCenter, pow(sphereRadius, 2));
//			cout << "Case 3!!" << endl;
			intersection(lineArc, s, back_inserter(intersections));
//			cout << "Case 3 ends!!" << endl;
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
		//cout << "Segment is of type 3!!" << endl;
		CGALPoint newPoint;
		
		float x = (A.x() + B.x()) / 2.0;
		float y = (A.y() + B.y()) / 2.0;
		float z = (A.z() + B.z()) / 2.0;
		newPoint = CGALPoint(x, y, z); 
		
		v = newPoint;
	}

	// update plc and DT
/*	cout << "DEBUG(Before plc update) !!" << endl;
	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
		cout << plc.point(pIter) << endl;
 	cout << "ENDS!!" << endl;
*/	updatePLCAndDT(v, missingSegmentHandle);
//	cout << "Going outside segment splitting function!!" << endl;
	return;
}


/*! \fn void CDTGenerator::recoverConstraintSegments()
    \brief Top-level routine for recovering constraint segments.

    Missing constriant segments of PLC are recovered by spliting the segments in such a way so that resulting sub-segments become _strongly Delaunay_.
*/
void CDTGenerator::recoverConstraintSegments()
{

	/*
	cout << "DEBUG(Before segment recovery)!!" << endl;
	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
		cout << plc.point(pIter) << endl;
*/
	DartHandle missingSegment;

	do
	{
		formMissingSegmentsQueue();
		size_t i = 0;
		while (missingSegmentQueue.size() != 0)
		{
			cout << "Segment recovery iteration: #" << i << endl;
			missingSegment = missingSegmentQueue.back();
			missingSegmentQueue.pop_back();
			splitMissingSegment(missingSegment);
			i++;
		}
	}while (missingSegmentQueue.size() != 0);

/*	////TEST
	cout << "Printing plc points after segment recovery: " << endl;
	for (LCC::One_dart_per_cell_range<2>::iterator fIter = plc.one_dart_per_cell<2>().begin(), fIterEnd = plc.one_dart_per_cell<2>().end(); fIter != fIterEnd; fIter++)
	{
		for (LCC::Dart_of_orbit_range<1>::iterator pIter = plc.darts_of_orbit<1>(fIter).begin(), pIterEnd = plc.darts_of_orbit<1>(fIter).end(); pIter != pIterEnd; pIter++)
			cout << plc.point(pIter) << " "; 
		cout << endl;
	}

	cout << "DEBUG(After segment recovery)!!" << endl;
	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
		cout << plc.point(pIter) << endl;
*/	
	return;
}


/*! \fn bool CDTGenerator::hasDegeneracyWithNeighbor(Delaunay::Cell_handle cHandle)
    \brief Tests whether input cell has any neighbor with which it forms cosphericality.
    \param [in] cHandle Handle to cell to be tested.
    \param [out] If there is a neighbor having degenerate condition with cHandle then it is returned here.
    \return True if cHandle has degenerate condition with a neighbor cell.    
*/
bool CDTGenerator::hasDegeneracyWithNeighbor(Delaunay::Cell_handle cHandle, size_t degenPartnerID)
{
	for (size_t neighborID = 0; neighborID < 4; neighborID++)
	{
		Delaunay::Cell_handle neighbor = cHandle->neighbor(neighborID);
		CGALPoint p[5];

		// initialize all other vertices of the 5-point set
		for (size_t i = 0, k = 0; i < 4; i++)
			p[k++] = cHandle->vertex(i)->point();
			
		p[4] = (neighbor->vertex(DT.mirror_index(cHandle, neighborID)))->point(); // the vertex opposite to cHandle in neighbor.
		if (side_of_bounded_sphere(p[0], p[1], p[2], p[3], p[4]) == ON_BOUNDARY)
		{
			degenPartnerID = neighborID;
			return true;
		}
		else
       			continue;
	}

	return false;       
}


/*! \fn CDTGenerator::correspondingVerticesInLCC(Delaunay::Cell_handle c1, size_t neighborID, DegenerateVertexSet &vertexSet)
 *  \brief Determines Dart handles to vertices involved in given degeneracy set.
 *  \param [in] c1 Handle to first cell.
 *  \param [in] neighborID Index of neighbor in c1.
 *  \param [out] vertexSet Vector of dart handles to vertices of LCC.
 */
void CDTGenerator::correspondingVerticesInLCC(Delaunay::Cell_handle c1, size_t neighborID, DegenerateVertexSet &vertexSet)
{
	// search for every vertex in cells inside plc 
	Delaunay::Cell_handle c2 = c1->neighbor(neighborID);
	size_t k = 0;
	for (size_t n = 0; n < 4; n++)
	{
		for (LCC::One_dart_per_cell_range<0>::iterator pHandle = plc.one_dart_per_cell<0>().begin(), pHandleEnd = plc.one_dart_per_cell<0>().end(); pHandle != pHandleEnd; pHandle++)
		{

			if (plc.point(pHandle) == c1->vertex(n)->point())
				vertexSet.vertHandles[k++] = pHandle;
			else
				continue;
		}
	}
	
	// add the vertex from c2 as well
	size_t neighborVertexIndex = DT.mirror_index(c1, neighborID);
	for (LCC::One_dart_per_cell_range<0>::iterator pHandle = plc.one_dart_per_cell<0>().begin(), pHandleEnd = plc.one_dart_per_cell<0>().end(); pHandle != pHandleEnd; pHandle++)
			if (plc.point(pHandle) == c2->vertex(neighborVertexIndex)->point())
				vertexSet.vertHandles[k++] = pHandle;
			else
				continue;
	
	
}


/*! \fn bool CDTGenerator::isVertexOnSegment(LCC::Dart_handle vertexHandle)
 *  \brief Tests whether input vertex is on a segment.
 *  \param [in] vertexHandle Dart handle to the input vertex.
 *  \return True if input vertex is on a segment. 
 */
bool CDTGenerator::isVertexOnSegment(LCC::Dart_handle vertexHandle)
{
	// take vertex 
	// check its neighborhood
	 
}


/*! \fn bool CDTGenerator::isVertexOnFacet(LCC::Dart_handle vertexHandle)
 *  \brief Tests whether input vertex is on a facet.
 *  \param [in] vertexHandle Dart handle to the input vertex.
 *  \return True if input vertex is on a facet. 
 */
bool CDTGenerator::isVertexOnFacet(LCC::Dart_handle vertexHandle)
{
}


/*! \fn bool CDTGenerator::isVertexPerturbable(LCC::Dart_handle perturbableVertexHandle)
 *  \brief Tests whether input vertex is _perturbable_.
 *  \param perturbableVertexHandle Handle to input vertex.
 *  \return True if input vertex is perturbable.
 */
bool CDTGenerator::isVertexPerturbable(LCC::Dart_handle perturbableVertexHandle)
{
	// Perturbable if symbolically perturbing it does not break consistency of PLC
	// Consistency:
		// Coplanarity
		// Collinearity
	//// STEPS:
	//////// 1. Parametric representation of problem.
	//////// 2. Check position of vertex: If on a segment or inside a facet then it is perturbable.
	//////// 3. If the vertex is an endpoint with more than one    	
	if (isVertexOnSegment(perturbableVertexHandle) || isVertexOnFacet(perturbableVertexHandle))
		return true;
	else // vertex at intersection of 3 or more facets, vertex on intersection of more than 3 segments
		return false;
	
	
}


/*! \fn bool CDTGenerator::hasPerturbableVertex(DegenerateVertexSet degenVertSet, LCC::Dart_handle &perturbableVertexHandle)
 *  \brief Tests whether given vertex set has atleast one perturbable vertex.
 *  \param [in] degenVertexSet Input vertex set of 5 degenerate vertices.
 *  \param [out] perturbableVertexHandle Handle to the perurbable vertex.
 *  \return True if there exists _atleast_ one degenerate vertex.
 */
bool CDTGenerator::hasPerturbableVertex(DegenerateVertexSet degenVertexSet, LCC::Dart_handle &perturbableVertexHandle)
{
	for (size_t i = 0; i < 5; i++) // test for each vertex
	{
		if (isVertexPerturbable(degenVertexSet.vertHandles[i]))
		{
			perturbableVertexHandle = degenVertexSet.vertHandles[i];
			return true;
		}
	}
	return false;
}


/*! \fn bool CDTGenerator::segmentSafePerturbable(LCC::Dart_handle vertexHandle)
 *  \brief Tests whether the perturbation of input vertex is _segment-safe_.
 *  \param [in] vertexHandle Handle to the vertex to be perturbed.
 *  \return True is input vertex is _segment-safe_ perturbable.
 */
bool CDTGenerator::segmentSafePerturbable(LCC::Dart_handle vertexHandle)
{
	
}


/*! \fn void CDTGenerator::perturbVertex(LCC::Dart_handle vertexToBePerturbed)
 *  \brief _Symbolically_ perturbs input vertex.
 *  \param [in] vertexToBePerturbed Handle to the vertex to be symbolically perturbed.
 */
void CDTGenerator::perturbVertex(LCC::Dart_handle vertexToBePerturbed)
{
	// express the vertex into 4d coordinates
	
	// add segment safe perturbation value to 4th coordinate of vertex
	for (LCC::One_dart_per_cell_range<1>::iterator segHandle = plc.one_dart_per_cell<1>().begin(), segHandleEnd = plc.one_dart_per_cell<1>().end(); segHandle != segHandleEnd; segHandle++)
	{
		// TODO: compute value of point perturbation.
	}
	// need to do optimization of perturbation value
	// compute Delaunay triangulation of perturbed vertex.
}


/*! \fn void CDTGenerator::removeLocalDegeneracies()
 *  \brief Breaks cosphericality condition(if exists) among any neighboring 5-points set in Delaunay triangulation. 
 */
void CDTGenerator::removeLocalDegeneracies()
{
	DegenerateVertexSet degenerateVertexSet;
	vector<DegenerateVertexSet> localDegeneracySet;
	size_t degenPartnerID;
	// Test for 5-point set of vertices which are cospherical
	for (Delaunay::Finite_cells_iterator cIter = DT.finite_cells_begin(), cIterEnd = DT.finite_cells_end(); cIter != cIterEnd; cIter++)
	{
		if (hasDegeneracyWithNeighbor(cIter, degenPartnerID))
		{
			correspondingVerticesInLCC(cIter, degenPartnerID, degenerateVertexSet); // handle to corresponding 5 vertices in PLC
			localDegeneracySet.push_back(degenerateVertexSet);		
		}
	}	
	
	// perturb a vertex in local degenerate vertex sets in 'segment safe' way.
	LCC::Dart_handle perturbableVertexHandle;
	for (vector<DegenerateVertexSet>::iterator dSetIter = localDegeneracySet.begin(), dSetIterEnd = localDegeneracySet.end(); dSetIter != dSetIterEnd; dSetIter++)
	{
		if (hasPerturbableVertex(*dSetIter, perturbableVertexHandle))
		{
			if (segmentSafePerturbable(perturbableVertexHandle))
				perturbVertex(perturbableVertexHandle); // symbolic perturbation
		}
	}
	// compute DT of modified vertices
	computeDelaunayTetrahedralization(-1); // TODO: Must consider the symbolic perturbation information
}


/*! \fn bool CDTGenerator::areFacetTetIntersecting(DartHandle tetHandle, DartHandle facetHandle)
 *  \brief Determines whether cell1 intersects cell2.
 *  \param [in] tetHandle Dart handle to first cell.
 *  \param [in] facetHandle Dart handle to second cell.
 *  \return true if tet anf facet intersect otherwise returns false.
 */
bool CDTGenerator::areFacetTetIntersecting(DartHandle tetHandle, DartHandle facetHandle)
{
	// represent facet as CGALTriangle and tet as CGALTeterahedron
	// call do_intersect(tet, tri)
	
	CGALTriangle tri;
	CGALTetrahedron tet;
	
	CGALPoint p1[3], p2[4];

	size_t i = 0;

	for (LCC::One_dart_per_incident_cell_range<0, 2>::iterator pIter = plc.one_dart_per_incident_cell<0, 2>(facetHandle).begin(), pIterEnd = plc.one_dart_per_incident_cell<0, 2>(facetHandle).end(); pIter != pIterEnd; pIter++)
		p1[i++] = plc.point(pIter);

	tri = CGALTriangle(p1[0], p1[1], p1[2]);

	i = 0;
	for (LCC::One_dart_per_incident_cell_range<0, 3>::iterator pIter = cdtMesh.one_dart_per_incident_cell<0, 3>(tetHandle).begin(), pIterEnd = cdtMesh.one_dart_per_incident_cell<0, 3>(tetHandle).end(); pIter != pIterEnd; pIter++)
		p2[i++] = plc.point(pIter);

	tet = CGALTetrahedron(p2[0], p2[1], p2[2], p2[3]);

	if (do_intersect(tet, tri))
		return true;
	else
		return false;
}

/*! \fn void CDTGenerator::computeMissingConstraintFacets(vector<DartHandle> &missingFacetList)
 *  \brief Computes list of constraint facets missing in current Delaunay triangulation.	
    \param [out] missingFacetList vector of Dart handles to missing constraint facets.	
 */
void CDTGenerator::computeMissingConstraintFacets(vector<DartHandle> &missingFacetList)
{
	// test which facets are not present in Delaunay triangulation 
	// add them to the missing facet list vector
	Delaunay::Vertex_handle v1, v2, v3;
	Delaunay::Cell_handle c;
	int i, j, k;

	for (LCC::One_dart_per_cell_range<2>::iterator fIter = plc.one_dart_per_cell<2>().begin(), fIterEnd = plc.one_dart_per_cell<2>().end(); fIter != fIterEnd; fIter++)
	{
		if (DT.is_vertex(plc.point(fIter), v1))
			if (DT.is_vertex(plc.point(plc.beta(fIter, 1)), v2))
				if (DT.is_vertex(plc.point(plc.beta(fIter, 1, 1)), v3))
					if (DT.is_facet(v1, v2, v3, c, i, j, k))
						continue;
		else
			missingFacetList.push_back(fIter); // facet is indeed missing from current output mesh.
	}
}


/*! \fn bool CDTGenerator::isNonStronglyDelaunayFacet(LCCWithDartInfo::Dart_handle d, LCCWithDartInfo lcc)
 *  \brief Tests whether input facet is not strongly Delaunay.
 *  \param [in] d Dart handle to the facet.
 *  \param [in] lcc Linear cell complex containing d.
 */
bool CDTGenerator::isNonStronglyDelaunayFacet(LCCWithDartInfo::Dart_handle d, LCCWithDartInfo& lcc)
{
	// a facet is non strongly Delaunay if there do not exist any circumsphere which does not include and other point on and inside it.
	// since local degeneracies are already removed we only need to check for input facet in Delaunay tetrahedralization of vetices of cavity.
	
	vector<CGALPoint> cavityPoints;
	for (LCCWithDartInfo::One_dart_per_cell_range<0>::iterator pIter = lcc.one_dart_per_cell<0>().begin(), pIterEnd = lcc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
		cavityPoints.push_back(lcc.point(pIter));

	Delaunay cavityDT;
	cavityDT.insert(cavityPoints.begin(), cavityPoints.end());
	
	Delaunay::Vertex_handle v1, v2, v3;
	Delaunay::Cell_handle c;
	int i, j, k;
	// test for facet
	if (cavityDT.is_vertex(lcc.point(d), v1))
		if (cavityDT.is_vertex(lcc.point(lcc.beta(d, 1)), v2))
			if (cavityDT.is_vertex(lcc.point(lcc.beta(d, 1, 1)), v3))
				if (cavityDT.is_facet(v1, v2, v3, c, i, j, k))
					return false;
	return true;
}


/*! \fn bool CDTGenerator::facetsHaveSameGeometry(LCC::Dart_handle fHandle, LCC lcc, LCCWithDartInfo::Dart_handle facetInCavity, LCCWithDartInfo cavityLCC)
 *   \brief Tests whether input facets from two different LCCs are _geometrically_ same.
 *  \param [in] fHandle Dart handle to the facet in LCC
 *  \param [in] lcc Linear cell complex constaining fHandle.
 *  \param [in] facetInCavity Dart handle to the facet in cavity.
 *  \param [in] cavityLCC LCC representation of cavity
 *  \return True if input facets have same geometry.   
 */
bool CDTGenerator::facetsHaveSameGeometry(LCC::Dart_handle fHandle, LCC& lcc, LCCWithDartInfo::Dart_handle facetInCavity, LCCWithDartInfo& cavityLCC)
{

	CGALPoint p[3];
	cout << "Inside facetsHaveSameGeometry!!" << endl;
	size_t i = 0;

	for (LCC::Dart_of_orbit_range<1>::iterator pHandleBegin = lcc.darts_of_orbit<1>(fHandle).begin(), pHandleEnd = lcc.darts_of_orbit<1>(fHandle).end(); pHandleBegin != pHandleEnd; pHandleBegin++)	
		p[i++] = lcc.point(pHandleBegin);

	LCC::Dart_handle d1 = lcc.make_triangle(p[0], p[1], p[2]);

	i = 0;
	for (LCCWithDartInfo::Dart_of_orbit_range<1>::iterator pIter = cavityLCC.darts_of_orbit<1>(facetInCavity).begin(),  pIterEnd = cavityLCC.darts_of_orbit<1>(facetInCavity).end(); pIter != pIterEnd; pIter++)
		p[i++] = cavityLCC.point(pIter);


	LCC::Dart_handle d2 = lcc.make_triangle(p[0], p[1], p[2]);
	if (lcc.are_facets_same_geometry(d1, d2))
		return true;
	else
		return false;
}


/*! \fn bool CDTGenerator::isFacetInCavity(LCC::Dart_handle fHandle, LCC lcc, LCCWithDartInfo::Dart_handle& correspondingFacetInCavity, LCCWithDartInfo cavityLCC)
 *  \brief Tests whether facet is in cavity
 *  \param [in] fHandle Dart handle to the facet in LCC
 *  \param [in] lcc Linear cell complex constaining fHandle.
 *  \param [out] correspondingFacetInCavity Returned facet in cavity.
 *  \param [in] cavityLCC LCC representation of cavity
 *  \return True if input facet is in cavity.
 */
bool CDTGenerator::isFacetInCavity(LCC::Dart_handle fHandle, LCC& lcc, LCCWithDartInfo::Dart_handle& correspondingFacetInCavity, LCCWithDartInfo& cavityLCC)
{
	cout << "Inside isFacetInCavity!!" << endl;
	size_t i = 0;
	for (LCCWithDartInfo::One_dart_per_cell_range<2>::iterator fIter = cavityLCC.one_dart_per_cell<2>().begin(), fIterEnd = cavityLCC.one_dart_per_cell<2>().end(); fIter != fIterEnd; fIter++)
	{
		if (facetsHaveSameGeometry(fHandle, lcc, fIter, cavityLCC))
		{
			correspondingFacetInCavity = fIter;
			return true;
		}
	}
	return false;
}


/*! \fn bool CDTGenerator::rayIntersectsFacet(CGALRay ray, LCCWithDartInfo::Dart_handle fHandle, LCCWithDartInfo& lcc)
 *  \brief Tests whether input ray intersect the facet in LCC.
 *  \param [in] ray Given ray.
 *  \param [in] fHandle Dart handle to facet in LCC.
 *  \param [in] lcc Linear cell complex containing fHandle.
 *  \return True if ray does intersect the given facet in LCC.
 */
bool CDTGenerator::rayIntersectsFacet(CGALRay ray, LCCWithDartInfo::Dart_handle fHandle, LCCWithDartInfo& lcc)
{
	// represent facet in triangle_3 form 
	CGALPoint p[3];
	size_t i = 0;
	for (LCCWithDartInfo::Dart_of_orbit_range<1>::iterator pIter = lcc.darts_of_orbit<1>(fHandle).begin(), pIterEnd = lcc.darts_of_orbit<1>(fHandle).end(); pIter != pIterEnd; pIter++)
		p[i++] = lcc.point(pIter);

	CGALTriangle triangle = CGALTriangle(p[0], p[1], p[2]);

	bool result = do_intersect(ray, triangle); 
	
	return result;
}	


/*! \fn bool CDTGenerator::isTetInsideCavity(Delaunay::Cell_handle ch, LCCWithDartInfo& cavityLCC)
 *  \brief Determines whether given 3-cell is inside Cavity boundary.
 *  \param [in] ch Cell handle pointing to the query tetrahedron.
 *  \param [in] cavityLCC LCC representation of cavity.
 *  \return True if given tetrahedron is inside cavityLCC.
 */
bool CDTGenerator::isTetInsideCavity(Delaunay::Cell_handle ch, LCCWithDartInfo& cavityLCC)
{

	CGALPoint circumcenter1 = circumcenter(ch->vertex(0)->point(), ch->vertex(1)->point(), ch->vertex(2)->point(), ch->vertex(3)->point());

	// generate a random ray from circumcenter
	Random r1 = Random(), r2 = Random(), r3 = Random(); // system time serves as seed.
	CGALPoint randomEndpoint(r1.get_bits<15>(), r2.get_bits<15>(), r3.get_bits<15>()); 
	CGALRay randomRay = CGALRay(circumcenter1, randomEndpoint);
	// find intersection with each face of cavity
	// if number of intersections are even then tetrahedron is outside
	// else, tetrahedron is outside.
	size_t nIntersections = 0;
	for (LCCWithDartInfo::One_dart_per_cell_range<2>::iterator fHandle = cavityLCC.one_dart_per_cell<2>().begin(), fHandleEnd = cavityLCC.one_dart_per_cell<2>().end(); fHandle != fHandleEnd; fHandle++)
	{
		if (rayIntersectsFacet(randomRay, fHandle, cavityLCC))
			nIntersections++;
		else 
			continue;
	}
	return (nIntersections % 2) ? false : true;
}


/*! \fn void CDTGenerator::recoverConstraintFacets()
 *  \brief Recovers constraint facets.
 */
void CDTGenerator::recoverConstraintFacets()
{
	// collect list of faces which are missing from DT
	// import DT to LCC
	// for each missing face:
	// 	compute 3-cells intersecting this facet
	// 	collect all these cells in a list
	//	compute cavity LCC:
	//		For each intersecting tetrahedron in cavity, extract the facets which should be part of cavityLCC
	//		Add those facets to cavityLCC
	// 	expand/verify cavity:     
	// 		Collect list of all facets L in cavityLCC which are not strongly Delaunay
	//		Remove that facet from L and cavityLCC
	//		Determine all other facets of the tetrahedra containing L:
	//			If that facet is already in L,  
	//				Yes, remove it from L and cavityLCC(if present)
	//				No, Add it to cavityLCC(it may be added in later iteration in L if it is not strongly Delaunay)
	//		Goto first step in cavity expansion	
	//	
	//	cavity retetrahedralization: 				
	//		Compute Delaunay tetrahedralization of vertices of cavityLCC
	//		For each tetrahedron in cavityLCC DT:
	//			Label it as inside/outside cavityLCC
	//		Sew all 'inside' tetrahedrons back to the hole creating in cdtMesh(LCC representation of output)

	cout << "Facet recovery starts..." << endl;	
	vector <DartHandle> missingConstraintFacets;
	vector <DartHandle> intersectingTets;
	LCCWithDartInfo cavityLCC;
	LCCWithDartInfo::Dart_handle cavityFaceHandle;

	computeMissingConstraintFacets(missingConstraintFacets); // list missing constraint facets
	DartHandle d = import_from_triangulation_3(cdtMesh, DT); // initialization of cdtMesh  
	
	// Remove infinite cells
	int infiniteVertexMark = cdtMesh.get_new_mark();
	markInfiniteVertexDart(d, cdtMesh, infiniteVertexMark);		
	vector <LCC::Dart_handle> cellsToBeRemoved;

	for (LCC::One_dart_per_cell_range<3>::iterator tetIter = cdtMesh.one_dart_per_cell<3>().begin(), tetIterEnd = cdtMesh.one_dart_per_cell<3>().end(); tetIter != tetIterEnd; tetIter++)
	{
		if (isInfinite(tetIter, cdtMesh, infiniteVertexMark, 3))
			cellsToBeRemoved.push_back(tetIter);
		else 
			continue;
	}
	cdtMesh.free_mark(infiniteVertexMark);

	for (vector<LCC::Dart_handle>::iterator cellIter = cellsToBeRemoved.begin(), cellIterEnd = cellsToBeRemoved.end(); cellIter != cellIterEnd; cellIter++)
		remove_cell<LCC, 3>(cdtMesh, *cellIter); // infinite cells removed from cdtMesh

	cout << "Infinite cells removed!!" << endl;

	while (missingConstraintFacets.size() != 0)
	{
		DartHandle missingFacetHandle = missingConstraintFacets.back(); // test for each missing facet
		missingConstraintFacets.pop_back();
		
		// compute cells intersecting THIS facet:
		intersectingTets.clear();
		for (LCC::One_dart_per_cell_range<3>::iterator cIter = cdtMesh.one_dart_per_cell<3>().begin(), cIterEnd = cdtMesh.one_dart_per_cell<3>().end(); cIter != cIterEnd; cIter++)
		{
			if (areFacetTetIntersecting(cIter, missingFacetHandle)) 
				intersectingTets.push_back(cIter);
		}
	
		int partOfIntersectingTetMark = cdtMesh.get_new_mark();
	
		if (partOfIntersectingTetMark == -1)
			exit(0);
	
		// Mark boundary facets of intersecting tets
		LCCWithDartInfo tempLCC; 
		LCCWithDartInfo::Dart_handle d;
		for (vector<LCC::Dart_handle>::iterator intersectingTetIter = intersectingTets.begin(); intersectingTetIterEnd = intersectingTets.end(); intersectingTetIter != intersectingTetIterEnd; intersectingTetIter++)
		{
			d = tempLCC.make_tetrahedron();
			// idetify the identical facets in the tet from cdtMesh and corresponding tet in tempLCC
			for (LCC::One_dart_per_incident_cell_range<2, 3>::iterator fIter1 = cdtMesh.one_dart_per_incident_cell<2, 3>(intersectingTetIter).begin(), fIterEnd1 = cdtMesh.one_dart_per_incident_cell<2, 3>(intersectingTetIter).end(); fIter1 != fIterEnd1; fIter1++)
				for (LCC::One_dart_per_incident_cell_range<2, 3>::iterator fIter2 = cdtMesh.one_dart_per_incident_cell<2, 3>(intersectingTetIter).begin(), fIterEnd2 = cdtMesh.one_dart_per_incident_cell<2, 3>(intersectingTetIter).end(); fIter2 != fIterEnd2; fIter2++)
 					if (areGeometricallySameFacets(fIter1, cdtMesh, fIter2, tempLCC)) 
						tempLCC.info<0>(fIter2) = fIter1;  // TODO: Is there any more efficient way that this to get dart to facet in cdtMesh and assign it to geometrically same facet in tempLCC??
		}
		tempLCC.sew3_same_facets();

		// mark faces which are at boundary
		int boundaryFacetMark = tempLCC.get_new_mark(); 
		if (boundaryFacetMark == -1)
		{
			cout << "BoundaryFacetMark: Free mark not available";
			exit(0);
		}	

		// copy the boundary facet to cavityLCC
		for (LCC::One_dart_per_incident_cell_range<2, 3>::iterator fIter = tempLCC.one_dart_per_incident_cell<2, 3>().begin(), fIterEnd = tempLCC.one_dart_per_incident_cell<2, 3>().end(); fIter != fIterEnd; fIterEnd++)
		{
			if (tempLCC.beta<3>(fIter) == tempLCC.null_dart_handle) // test for getting boundary facet
			{
				d = cavityLCC.make_triangle(tempLCC.point(fIter), tempLCC.point(tempLCC.beta<1>(fIter), tempLCC.point(tempLCC.beta<1, 1>(fIter));
				cavityLCC.info(d) = tempLCC.info(fIter); // stores handle to the facet in original mesh	
			}
		}
		
		/*	for (vector<DartHandle>::iterator intersectingTetIter = intersectingTets.begin(), intersectingTetIterEnd = intersectingTets.end(); intersectingTetIter != intersectingTetIterEnd; intersectingTetIter++)
		{
			for (LCC::One_dart_per_incident_cell_range<2, 3>::iterator faceIter = cdtMesh.one_dart_per_incident_cell<2, 3>(*intersectingTetIter).begin(), faceIterEnd = cdtMesh.one_dart_per_incident_cell<2, 3>(*intersectingTetIter).end(); faceIter != faceIterEnd; faceIter++)
			{
				cdtMesh.mark(faceIter, partOfIntersectingTetMark);
				if (cdtMesh.beta<3>(faceIter) != cdtMesh.null_dart_handle) // since a face is shared by 2 3-cells
					cdtMesh.mark(cdtMesh.beta<3>(faceIter), partOfIntersectingTetMark);
			}
		}
	
		//// add marked faces to cavityLCC
		for (LCC::Dart_range::iterator fHandle = cdtMesh.darts().begin(), fHandleEnd = cdtMesh.darts().end(); fHandle != fHandleEnd; fHandle++)
			if (cdtMesh.is_marked(fHandle, partOfIntersectingTetMark))
			{
				CGALPoint p[3];

				size_t i = 0;
				for (LCC::Dart_of_orbit_range<1>::iterator pIter = cdtMesh.darts_of_orbit<1>(fHandle).begin(), pIterEnd = cdtMesh.darts_of_orbit<1>(fHandle).end(); pIter != pIterEnd; pIter++)
					p[i++] = cdtMesh.point(pIter);
				
				cavityFaceHandle = cavityLCC.make_triangle(p[0], p[1], p[2]);
				for (LCCWithDartInfo::Dart_of_orbit_range<1>::iterator pIter = cavityLCC.darts_of_orbit<1>(cavityFaceHandle).begin(), pIterEnd = cavityLCC.darts_of_orbit<1>(cavityFaceHandle).end(); pIter != pIterEnd; pIter++)	
					cavityLCC.info<0>(pIter) = fHandle;// handle to that facet in original mesh 				
			}
		cdtMesh.free_mark(partOfIntersectingTetMark); 
*/
		// Add only those faces which are on boundary
		 
		 
		//// remove intersecting tets from cdtMesh
		for (vector<LCC::Dart_handle>::iterator tetIter = intersectingTets.begin(), tetIterEnd = intersectingTets.end(); tetIter != tetIterEnd; tetIter++)
			remove_cell<LCC, 3>(cdtMesh, *tetIter); // TODO: shouldn't it be at the end?
		cout << "Intersecting tets removed from mesh!!" << endl;

		//// sew 2-cells at boundaries
		sew2CellsWithDartInfoFromEdge(cavityLCC); // cavity is created
		
		cout << "Cavity verification started!!" << endl;
		// CAVITY VERIFICATION/EXPANSION:
		//// create queue of non strongly Delaunay faces in cavityLCC
		vector<LCCWithDartInfo::Dart_handle> nonStronglyDelaunayFacetsInCavity;
		LCCWithDartInfo::Dart_handle correspondingFacetInCavity;

		do
		{
			//// initialize vector
			for (LCCWithDartInfo::One_dart_per_cell_range<2>::iterator nonStrongFaceIter = cavityLCC.one_dart_per_cell<2>().begin(), nonStrongFaceIterEnd = cavityLCC.one_dart_per_cell<2>().end(); nonStrongFaceIter != nonStrongFaceIterEnd; nonStrongFaceIter++)	
				if (isNonStronglyDelaunayFacet(nonStrongFaceIter, cavityLCC))  
					nonStronglyDelaunayFacetsInCavity.push_back(nonStrongFaceIter); // handle to facets in cavityLCC
				else 
					continue;	
			cout << "Found non-strongly Delaunay facets in cavity!!" << endl;
			//// cavity expansion
			while (nonStronglyDelaunayFacetsInCavity.size() != 0)	
			{
				LCCWithDartInfo::Dart_handle nonStronglyDelaunayFace = nonStronglyDelaunayFacetsInCavity.back();
				nonStronglyDelaunayFacetsInCavity.pop_back();
					
				LCC::Dart_handle exteriorCellSharingNonDelaunayFacet = cavityLCC.info<0>(nonStronglyDelaunayFace);  
				//// Explore all faces of this cell
				size_t i = 0;
				for (LCC::One_dart_per_incident_cell_range<2, 3>::iterator facetInCellHandle = cdtMesh.one_dart_per_incident_cell<2, 3>(exteriorCellSharingNonDelaunayFacet).begin(), facetInCellEndHandle = cdtMesh.one_dart_per_incident_cell<2, 3>(exteriorCellSharingNonDelaunayFacet).end(); facetInCellHandle != facetInCellEndHandle; facetInCellHandle++) 
				{	
					cout << "Iteration: " << i++ << endl;
					size_t facetLocation;
					cout << "Before isFacetInCavity!!" << endl;
					if (isFacetInCavity(facetInCellHandle, cdtMesh, correspondingFacetInCavity, cavityLCC)) 
					{
						if (correspondingFacetInCavity != cavityLCC.null_dart_handle)
							remove_cell<LCCWithDartInfo, 2>(cavityLCC, correspondingFacetInCavity);
						else
						{
							cout << "Dart returned by isFacetInCavity is NULL!!";
							exit(0);
						}
					
					}	
					else
					{
						//cout << "In else!!" << endl;
						CGALPoint p[3];
						size_t i = 0;
						for (LCC::Dart_of_orbit_range<1>::iterator pIter = cdtMesh.darts_of_orbit<1>(facetInCellEndHandle).begin(), pIterEnd = cdtMesh.darts_of_orbit<1>(facetInCellEndHandle).end(); pIter != pIterEnd; pIter++)
							p[i++] = cdtMesh.point(pIter);
				
						cavityLCC.make_triangle(p[0], p[1], p[2]);
						sew2CellsWithDartInfoFromEdge(cavityLCC);
						//cout << "After sewing!!" << endl;
					}
				}
			}
			cout << "Cavity expanded, if required!!" << endl;
			
		}while (nonStronglyDelaunayFacetsInCavity.size() != 0);
		cout << "Cavity verification complete!!" << endl;
		cout << "Cavity retetrahedralization started!!" << endl;
		// CAVITY RETETRAHEDRALIZATION		
		vector<pair<CGALPoint, LCC::Dart_handle> > cavityVertices;
		
		for (LCCWithDartInfo::One_dart_per_cell_range<0, 3>::iterator pIter = cavityLCC.one_dart_per_cell<0, 3>().begin(), pIterEnd = cavityLCC.one_dart_per_cell<0, 3>().end(); pIter != pIterEnd; pIter++)	
			cavityVertices.push_back(make_pair(cavityLCC.point(pIter), cavityLCC.info<0>(pIter)));
	
		Delaunay cavityDelaunay;
		cavityDelaunay.insert(cavityVertices.begin(), cavityVertices.end()); // DT of cavity vertices
	
		// mark each cell of DT as either inside/outside cavity
		for (Delaunay::Finite_cells_iterator cIter = cavityDelaunay.finite_cells_begin(), cIterEnd = cavityDelaunay.finite_cells_end(); cIter != cIterEnd; cIter++)
		{
			if (isTetInsideCavity(cIter, cavityLCC))
			{
				// create identical 3-cell in cdtMesh
				CGALPoint p[4];
				for (size_t i = 0; i < 4; i++)
					p[i] = ((*cIter).vertex(i))->point();
				cdtMesh.make_tetrahedron(p[0], p[1], p[2], p[3]); 
			}
			else
				continue;
		}
			
		cdtMesh.sew3_same_facets();
		cout << "Cavity retetrahedralization complete!!" << endl;
	}

	cout << "Constraint facets recovered!!" << endl;

}


/*! \fn bool CDTGenerator::rayIntersectsPLCFacet(CGALRay randomRay, LCC::Dart_handle fHandle)
 *  \brief Tests whether input ray intersects the PLC facet. 
 *  \param [in] randomRay Random ray 
 *  \param [in] fHandle Dart handle to the facet. 
 *  \return True if input ray intersects the facet of PLC.
 */
void CDTGenerator::countRayPLCFacetIntersections(CGALRay randomRay, LCC::Dart_handle fHandle, size_t &nIntersections)
{
	// represent facet in triangle_3 form 
	CGALPoint p[4];
	size_t i = 0;

	for (LCC::One_dart_per_incident_cell_range<0, 2>::iterator pIter = plc.one_dart_per_incident_cell<0, 2>(fHandle).begin(), pIterEnd = plc.one_dart_per_incident_cell<0, 2>(fHandle).end(); pIter != pIterEnd; pIter++)
	{
		p[i++] = plc.point(pIter);
//		cout << "Point: " << plc.point(pIter) << endl;
		cout << "Point: " << p[i - 1] << endl;
	}
	CGALTriangle triangle;
	bool result;

	triangle = CGALTriangle(p[0], p[1], p[2]);
	if (triangle.is_degenerate())
	{
		cout << "Degenerate triangle 1...discarding!!" << endl;
	}
	else
	{
		if (do_intersect(randomRay, triangle)) 
		{
			cout << "Intersects!!" << endl; 
			nIntersections++;
		}
	}
	if (i == 4) // facet has 2 triangles
	{
		triangle = CGALTriangle(p[0], p[2], p[3]);
		if (triangle.is_degenerate())
		{
			cout << "Degenerate triangle 2...discarding!!" << endl;
		}
		else
		{
			if (do_intersect(randomRay, triangle))
			{
				nIntersections++;			
				cout << "Intersects!!" << endl;
			}
		}
	}
	//	cout << "Source point of the ray is: " << randomRay.source() << endl;

}


/*! \fn bool CDTGenerator::isCellOutsidePLC(LCC::Dart_handle cellHandle)
 *  \brief Tests whether input cell is outside PLC
 *  \param [in] cellHandle Dart handle for the input cell.
 *  \return True if the input cell is outside PLC.
 */
bool CDTGenerator::isCellOutsidePLC(LCC::Dart_handle cellHandle)
{
	// endpoint 1
	size_t i = 0;
	CGALPoint p[4];
	for (LCC::One_dart_per_incident_cell_range<0, 3>::iterator pIter = cdtMesh.one_dart_per_incident_cell<0, 3>(cellHandle).begin(), pIterEnd = cdtMesh.one_dart_per_incident_cell<0, 3>(cellHandle).end(); pIter != pIterEnd; pIter++)
		p[i++] = cdtMesh.point(pIter); 

	CGALPoint testTetCentroid = centroid(p[0], p[1], p[2], p[3]);

	CGALPoint plcCentroid;	
	float x = 0.0, y = 0.0, z = 0.0;
	size_t nPoints = 0;

	// endpoint 2
	for (LCC::One_dart_per_cell_range<0>::iterator pIter = plc.one_dart_per_cell<0>().begin(), pIterEnd = plc.one_dart_per_cell<0>().end(); pIter != pIterEnd; pIter++)
	{
		x += plc.point(pIter).x();
		y += plc.point(pIter).y();
		z += plc.point(pIter).z();
		nPoints++;
	}
	plcCentroid = CGALPoint(x / nPoints, y / nPoints, z / nPoints);

	// compute random ray
	CGALRay randomRay = CGALRay(testTetCentroid, plcCentroid);
	 
	if (randomRay.is_degenerate())
	{
		cout << "Generated ray is degenerate!!";
		exit(0);
	}
	cout << "Centroid of plc is:  " << plcCentroid << endl;
	cout << "Centroid of test tet is:  " << testTetCentroid << endl;

	size_t nIntersections = 0;
	for (LCC::One_dart_per_cell_range<2>::iterator fHandle = plc.one_dart_per_cell<2>().begin(), fHandleEnd = plc.one_dart_per_cell<2>().end(); fHandle != fHandleEnd; fHandle++)
		//if (cdtMesh.beta<2>(fHandle) == NULL) // boundary facet
		countRayPLCFacetIntersections(randomRay, fHandle, nIntersections);
		
	cout << "Total number of intersections: " << nIntersections << endl;
	return (nIntersections % 2 == 0) ? true : false;
}



/*! \fn void CDTGenerator::removeExteriorTetrahedrons()
 *  \brief Removes tetrahedrons from cdtMesh which are exterior wrt. input PLC.
 */
void CDTGenerator::removeExteriorTetrahedrons()
{
	// Ray-LCC intersection approach, if even number of intersections, then outside else inside. 
	vector<LCC::Dart_handle> exteriorCellsList;

	//// For each 3-cell
	for (LCC::One_dart_per_cell_range<3>::iterator cellIter = cdtMesh.one_dart_per_cell<3>().begin(), cellIterEnd = cdtMesh.one_dart_per_cell<3>().end(); cellIter != cellIterEnd; cellIter++)
	{
		//// Test if is it exterior cell
		if (isCellOutsidePLC(cellIter))
			exteriorCellsList.push_back(cellIter);
	}
	
	size_t nCells = 0;
	for (LCC::One_dart_per_cell_range<3>::iterator cellIter = cdtMesh.one_dart_per_cell<3>().begin(), cellIterEnd = cdtMesh.one_dart_per_cell<3>().end(); cellIter != cellIterEnd; cellIter++)
		nCells++;
	
	cout << "Number of cells before removal: " << nCells << endl;
	
	//// remove cells marked as exterior
	for (vector<LCC::Dart_handle>::iterator cellsToBeRemovedIter = exteriorCellsList.begin(), cellsToBeRemovedIterEnd = exteriorCellsList.end(); cellsToBeRemovedIter != cellsToBeRemovedIterEnd; cellsToBeRemovedIter++)
		remove_cell<LCC, 3>(cdtMesh, *cellsToBeRemovedIter);			

	cout << "Number of cells to be removed: " << exteriorCellsList.size() << endl;
	//nCells = 0;
	//for (LCC::One_dart_per_cell_range<3>::iterator cellIter = cdtMesh.one_dart_per_cell<3>().begin(), cellIterEnd = cdtMesh.one_dart_per_cell<3>().end(); cellIter != cellIterEnd; cellIter++)
	//	nCells++;
	
	//cout << "Number of cells after removal: " << nCells << endl;
	// write mesh to PLY file
	//writePLYOutput(NULL, cdtMesh, "../../data/outputMesh.ply");
}


/*! \fn void CDTGenerator::generate()
    \brief Public interface for CDTGenerator class.
 */
void CDTGenerator::generate()
{
	readPLCInput();
	computeDelaunayTetrahedralization(-1);
	recoverConstraintSegments();
//	removeLocalDegeneracies();
	recoverConstraintFacets();
	removeExteriorTetrahedrons(); // removes tetrahedrons from cdtMesh which are outside input PLC

}



