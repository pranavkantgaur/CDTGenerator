#include <iostream>
#include "rply/rply.h"
#include <stdlib.h>

using namespace std;


int main()
{

	float p[5][3];

	for (unsigned int k = 0; k < 5; k++)
	{
		p[k][0] = 1.0;
		p[k][1] = 2.0;
		p[k][2] = 3.0;
	}

	p_ply pointArray;

	if((pointArray = ply_create("output.ply", PLY_ASCII, NULL, 0, NULL)) == NULL)
	{
		cout << "\nCannot write Mesh output!!";
		exit(0);
	}
	 
	else
	{
		ply_add_element(pointArray, "vertex", 5); 
	      	ply_add_scalar_property(pointArray, "x", PLY_FLOAT);
		ply_add_scalar_property(pointArray, "y", PLY_FLOAT);
		ply_add_scalar_property(pointArray, "z", PLY_FLOAT);
	
		if (!ply_write_header(pointArray))
			cout << "\nHeader not writen" << flush;

		for (unsigned int n = 0; n < 5; n++)
			for (unsigned int m = 0; m < 3; m++)
				ply_write(pointArray, p[n][m]);

	
		ply_close(pointArray);
	}

	return 0;
}


