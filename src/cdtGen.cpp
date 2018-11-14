#include "cdtGen.h"

static size_t pointId = 0;
static float tempPoint[3];
static size_t dimensionId = 0;
vector<CGALPoint> plcVertexVector;
vector<TriangleWithIndices> plcFaceVector;


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
        static TriangleWithIndices tempFace;
	
	ply_get_argument_property(argument, NULL, &length, &value_index);

        switch (value_index) 
	{
        	case 0:
	        case 1: 
        		tempFace.pointIDs[pointId++] = ply_get_argument_value(argument);
		        break;
        	case 2:	
			tempFace.pointIDs[pointId] = ply_get_argument_value(argument);
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
bool CDTGenerator::areGeometricallySameSegments(DartHandle& d1, DartHandle& d2, LCC& lcc)
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






/*! \fn void CDTGenerator::readPLCInput()
    \brief Reads PLC from input(only PLY supported currently) file.

    Reads vertex and face information(from respective vectors) and initializes corresponding linear cell complex _plc_.
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
			vertexIds[k] = plcFaceVector[n].pointIDs[k];
	
		for (size_t i = 0; i < 3; i++)
			trianglePoints[i] = CGALPoint(plcVertexVector[vertexIds[i]].x(), plcVertexVector[vertexIds[i]].y(), plcVertexVector[vertexIds[i]].z());
			
		plc.make_triangle(trianglePoints[0], trianglePoints[1], trianglePoints[2]);
	}
	
	// sew facets sharing edge
	sew2CellsFromEdge(plc);
}



/*! \fn void CDTGenerator::generate()
    \brief Public interface for CDTGenerator class.
 */
void CDTGenerator::generate()
{
	readPLCInput();
//	computeDelaunayTetrahedralization(-1);
//	recoverConstraintSegments();
//	removeLocalDegeneracies();
/*	cout << "Skipping explicit local degeneracy removal, CGAL performs symbolic perturbation by default!!" << endl;
	recoverConstraintFacets();
	removeExteriorTetrahedrons(); // removes tetrahedrons from cdtMesh which are outside input PLC
*/
	return;
}

