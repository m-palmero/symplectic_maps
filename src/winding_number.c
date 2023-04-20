// --- Function Time Series ---
// Evolve a i.c. and save the time series evolution of a specific coodinate
# include "dynamics.h"

extern char filename1[T], filename2[T];
extern int iter; 
extern bool standard_map;
extern bool boozer_map;
extern bool ullmann_map;

int compute_winding_number(MAP *m)
{	
	FILE *out1 = fopen(filename1, "w");
	FILE *out2 = fopen(filename2, "w");

	double x0[2], x1[2];
	double y, x, w;
	double deltay;
	double y_max;

    x0[0] = m->x0_center;
    x0[1] = m->y0_center;

	//y_max = 4.5;
	//ny = 2000;

	//deltay = fabs(y_max - m->y0_center) / ((double)ny - 1);

	x = m->x0_center;
	y = m->y0_center;
	//for (int j = 0; j < ny; j++)
	//{
		x0[0] = x;
    	x0[1] = y;
		fprintf(out2, "%1.16e %1.16e\n", x0[0],x0[1]);
		for (int k = 0; k < iter; k++)
    	{	
        	if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
			if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
			if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
			x0[0] = x1[0];
			x0[1] = x1[1];
			//fprintf(out1, "%1.16e %1.16e\n", x1[0],x1[1]);
		}
		//if (ny > 1) y += deltay;
		w = fabs(x1[0] - m->x0_center)/(iter * 2 * M_PI);
		printf("winding number = %f\n", w);
		//fprintf(out1, "%1.16e %1.16e\n", y, w);
	//}

	fclose(out1);
	fclose(out2);
    
    
    return 0;
}