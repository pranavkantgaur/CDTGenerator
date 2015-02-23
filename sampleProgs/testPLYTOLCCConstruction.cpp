/*
 * This program reads a PLY file containing description of a polyhedral domain and generates corresponding LCC representation
 */

void readPLYDataset()
{
	initialize POLYGONS;
}

void generateLCC()
{
	for (each polygon p in POLYGONS)
		lcc.make_triangle(p);
	
	for (triangle t in lcc)
	{
		for (dart handle d1 of t)
		{
			if (d1.mark == false) // not sewed till now
				for (triangle t' in lcc && t' != t)
				{
					for (dart handle d2 of t')
						if (shareSegment(d1,d2))
						{
							sew<2>(d1, d2);
							d1.mark = d2.mark = true; // mark them as sewed
						}
				}
		}
	}
}

int main()
{
	readPLYDataset();
	generateLCC();
	return 0;
}
