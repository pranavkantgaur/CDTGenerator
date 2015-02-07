/* 
 * Minimal test for testing order of vertices returned by finite_facet_iterator
 */
#include <vector>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>


using namespace std;
using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;
typedef Triangulation_vertex_base_with_info_3<unsigned int, K> vbWithInfo;
typedef Triangulation_data_structure_3<vbWithInfo> tdsWithInfo;
typedef Delaunay_triangulation_3<K, tdsWithInfo, Fast_location> Delaunay;
typedef Delaunay::Point Point;


int main()
{
	vector<pair<Point, unsigned int> > samplePoints;

	samplePoints.push_back(pair<Point, unsigned int>(Point(0, 0, 0), 0)); // points of a cube
	samplePoints.push_back(pair<Point, unsigned int>(Point(0, 0, 1), 1));
	samplePoints.push_back(pair<Point, unsigned int>(Point(0, 1, 0), 2));
	samplePoints.push_back(pair<Point, unsigned int>(Point(0, 1, 1), 3));
	samplePoints.push_back(pair<Point, unsigned int>(Point(1, 0, 0), 4));
	samplePoints.push_back(pair<Point, unsigned int>(Point(1, 0, 1), 5));
	samplePoints.push_back(pair<Point, unsigned int>(Point(1, 1, 0), 6));
	samplePoints.push_back(pair<Point, unsigned int>(Point(1, 1, 1), 7));

	Delaunay DT(samplePoints.begin(), samplePoints.end());

	for (Delaunay::Finite_facets_iterator fIter = DT.finite_facets_begin(); fIter != DT.finite_facets_end(); fIter++)
	{
		for (unsigned int i = 0; i < 4; i++)
			if (i != fIter->second)
			{
				cout << fIter->first->vertex(i)->info() << "\t" << endl;
			}
		cout << "\n";
	}

	return 0;
}

