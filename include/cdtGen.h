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

#include "rply.h"


#define Pi 22.0/7.0
#define INVALID_VALUE -1.0f // used in context of distances 
using namespace std;
using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;
typedef Linear_cell_complex_traits<3, K> Traits;
typedef Linear_cell_complex<3, 3, Traits> LCC;
typedef LCC::Dart_handle DartHandle;

struct MyDartInfo
{
	template<class Refs>
	struct Dart_wrapper
	{
		typedef CGAL::Dart<3, Refs > Dart;
		typedef Cell_attribute_with_point<Refs, DartHandle > VertexAttribute;
		typedef cpp11::tuple<VertexAttribute> Attributes;
	};	
};

struct MyIntInfo
{
	template<class Refs>
	struct Dart_wrapper
	{
		typedef CGAL::Dart<3, Refs > Dart;
		typedef Cell_attribute_with_point<Refs, size_t > VertexAttribute;
		typedef cpp11::tuple<VertexAttribute> Attributes;
	};	
};

typedef Linear_cell_complex<3, 3, Traits, MyDartInfo> LCCWithDartInfo;
typedef Linear_cell_complex<3, 3, Traits, MyIntInfo> LCCWithIntInfo;
typedef LCCWithDartInfo::Dart_handle DartHandleWithDartInfo;
typedef Triangulation_vertex_base_with_info_3<DartHandle, K> Vb; 
typedef Triangulation_data_structure_3<Vb> Tds;
typedef Delaunay_triangulation_3<K, Tds, Fast_location> Delaunay;
typedef Delaunay::Vertex_handle VertexHandle;
typedef Delaunay::Cell_handle CellHandle;
typedef Point_3<K> CGALPoint;
typedef Circle_3<K> CGALCircle;
typedef Sphere_3<K> CGALSphere;
typedef Tetrahedron_3<K> CGALTetrahedron;
typedef Triangle_3<K> CGALTriangle;
typedef Ray_3<K> CGALRay;
typedef Segment_3<K> CGALSegment;
typedef Exact_spherical_kernel_3 SK;
typedef Line_arc_3<SK> CGALSphericalLineArc;
typedef Sphere_3<SK> CGALSphericalSphere;
typedef Segment_3<SK> CGALSphericalSegment; 
typedef Point_3<SK> CGALSphericalPoint;

/*! \class Triangle
    \brief Represents triangle 
  
    Triangle class stores indices of vertices. 
*/
class Triangle
{
	public:
		size_t pointIds[3]; /*< Indices of points */
};


class CDTGenerator
{
	public:
		void generate();
	protected:

		LCC plc; /*!< Input piecewise linear cell complex representing the input */
		Delaunay DT; /*!< Intermidiate structure used for storing Delaunay tetrahedralization */
		LCC cdtMesh; /*!< Output mesh */
		
		void markInfiniteVertexDart(LCCWithIntInfo::Dart_handle, LCCWithIntInfo&, int);
		bool isInfinite(LCCWithIntInfo::Dart_handle, LCCWithIntInfo, int, unsigned int);
		void recoverConstraintSegments();
		void splitMissingSegment(DartHandle);
		void updatePLCAndDT(CGALPoint&, DartHandle);
		float computeSegmentLength(CGALPoint&, CGALPoint&);
		unsigned int determineSegmentType(DartHandle);
		bool isVertexAcute(DartHandle);
		float computeAngleBetweenSegments(DartHandle, DartHandle);
		float vectorMagnitude(DartHandle);
		float dotProduct(DartHandle, DartHandle);
		void computeReferencePoint(CGALPoint*, DartHandle);
		unsigned int computeCircumradius(CGALPoint&, CGALPoint&, CGALPoint&);
		void formMissingSegmentsQueue(vector<DartHandle>&);
		void computeDelaunayTetrahedralization();
		void writePLYOutput(LCCWithIntInfo::Dart_handle, LCCWithIntInfo&, string); 
		void readPLCInput();
		bool areGeometricallySameSegments(DartHandle, DartHandle);
};
