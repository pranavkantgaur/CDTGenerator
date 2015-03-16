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
typedef Ray_3<K> CGALRay;

typedef Exact_spherical_kernel_3 SK;
typedef Line_arc_3<SK> CGALSphericalSegment;
typedef Sphere_3<SK> CGALSphericalSphere;
typedef Segment_3<SK> CGALSphericalSegment; 
typedef Point_3<SK> CGALSphericalPoint;
