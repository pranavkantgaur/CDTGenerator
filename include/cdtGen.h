#include <iostream>
#include <vector>
#include <math.h>
#include <map>
#include <unordered_set>

#include <CGAL/Kernel/global_functions.h>
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
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>
#include <utility>


#include "rply.h"


#define Pi 22.0/7.0
#define INVALID_VALUE -1.0f // used in context of distances 
using namespace std;
using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;
typedef Linear_cell_complex_traits<3, K> Traits;
typedef Linear_cell_complex<3, 3, Traits> LCC;
typedef LCC::Dart_handle DartHandle;
typedef boost::tuple<Point_3,int>                           Point_and_int;
typedef CGAL::Random_points_in_cube_3<Point_3>              Random_points_iterator;
typedef CGAL::Search_traits_3<Kernel>                       Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int,
  CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
  Traits_base>                                              Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>          K_neighbor_search;
typedef K_neighbor_search::Tree                             Tree;
typedef K_neighbor_search::Distance                         Distance;

/*
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
*/

struct MyDartInfo
{
	template<class CMap>
	struct Dart_wrapper
	{
		typedef CGAL::Dart<3, CMap > Dart;
		typedef Cell_attribute<CMap, DartHandle> Facet_attribute;
		typedef Cell_attribute_with_point<CMap> Vertex_attribute;
		typedef	cpp11::tuple<Vertex_attribute, void, Facet_attribute> Attributes;
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
typedef std::list<CGALTriangle>::iterator Iterator;
typedef AABB_triangle_primitive<K, Iterator> Primitive;
typedef AABB_traits<K, Primitive> AABB_triangle_traits;
typedef AABB_tree<AABB_triangle_traits> aabbTree;

/*! \class Triangle
    \brief Represents triangle 
  
    Triangle class stores indices of vertices. 
*/
class Triangle
{
	public:
		size_t pointIds[3]; /*< Indices of points */
};


class DegenerateVertexSet
{
	public:
		LCC::Dart_handle vertHandles[5];
};


class CDTGenerator
{
	public:
		void generate();
	protected:

		LCC plc; /*!< Input piecewise linear cell complex representing the input */
		Delaunay DT; /*!< Intermidiate structure used for storing Delaunay tetrahedralization */
		LCC cdtMesh; /*!< Output mesh */
		vector<LCC::Dart_handle> missingSegmentQueue;

		void markInfiniteVertexDart(LCCWithIntInfo::Dart_handle, LCCWithIntInfo&, int);
		void markInfiniteVertexDart(LCC::Dart_handle, LCC&, int);

		bool isInfinite(LCCWithIntInfo::Dart_handle&, const LCCWithIntInfo&, int, size_t);
		bool isInfinite(LCC::Dart_handle&, const LCC&, int, size_t);
		void recoverConstraintSegments();
		void splitMissingSegment(DartHandle&);
		void updatePLCAndDT(CGALPoint&, DartHandle&);
		float computeSegmentLength(CGALPoint&, CGALPoint&);
		size_t determineSegmentType(DartHandle&);
		bool isVertexAcute(DartHandle);
		float dotProduct(DartHandle&, DartHandle&);
		void computeReferencePoint(CGALPoint*, DartHandle&);
		size_t computeCircumradius(CGALPoint&, CGALPoint&, CGALPoint&);
		void formMissingSegmentsQueue();//vector<DartHandle>&);
		void computeDelaunayTetrahedralization(int);
		void copyLCCToLCCWithIntInfo(LCC &, LCCWithIntInfo &);
		void writePLYOutput(LCCWithIntInfo::Dart_handle, LCCWithIntInfo&, string); 
		void writePLYOutput(LCC::Dart_handle, LCC&, string); 

		void readPLCInput();

		bool areGeometricallySameSegments(DartHandle&, DartHandle&, LCC&);
		bool areFacetsGeometricallySame(LCC::Dart_handle&, LCC&, LCC::Dart_handle&, LCC&);
		bool areGeometricallySameSegmentsWithDartInfo(LCCWithDartInfo::Dart_handle&, LCCWithDartInfo::Dart_handle&, LCCWithDartInfo&);
		void sew2CellsFromEdge(LCC &);
		void sew2CellsWithDartInfoFromEdge(LCCWithDartInfo &);

		// local degeneracy removal
		bool hasDegeneracyWithNeighbor(Delaunay::Cell_handle, size_t); 
		void correspondingVerticesInLCC(Delaunay::Cell_handle, size_t, DegenerateVertexSet&); 
		bool isVertexOnSegment(LCC::Dart_handle);
		bool isVertexOnFacet(LCC::Dart_handle);
		bool isVertexPerturbable(LCC::Dart_handle);
		bool hasPerturbableVertex(DegenerateVertexSet, LCC::Dart_handle&);
		bool segmentSafePerturbable(LCC::Dart_handle);
		void perturbVertex(LCC::Dart_handle);
		void removeLocalDegeneracies();

		// constraint facet recovery
		bool areFacetTetIntersecting(DartHandle&, DartHandle&); // first 2 arguments specify the dimension of first and second cells respectively.
		void computeMissingConstraintFacets(vector<DartHandle>&);
		void recoverConstraintFacets();
		bool isNonStronglyDelaunayFacet(LCCWithDartInfo::Dart_handle&, LCCWithDartInfo&); //TODO
		bool facetsHaveSameGeometry(LCC::Dart_handle&, LCC&, LCCWithDartInfo::Dart_handle&, LCCWithDartInfo&); // TODO
		bool isFacetInCavity(LCC::Dart_handle&, LCC&, LCCWithDartInfo::Dart_handle&, LCCWithDartInfo&);
	
		bool isTetInsideCavity(Delaunay::Cell_handle, LCCWithDartInfo&);
		bool rayIntersectsFacet(CGALRay&, LCCWithDartInfo::Dart_handle&, LCCWithDartInfo&);
		void countRayPLCFacetIntersections(CGALRay&, LCC::Dart_handle&, size_t &);
		bool isCellOutsidePLC(LCC::Dart_handle&, CGALPoint&);
		void removeExteriorTetrahedrons();

};
