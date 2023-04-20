# include "dynamics.h"

extern char filename1[T];
extern double parameter;
extern int iter;
extern int nx;
extern int ny;
extern bool standard_map;
extern bool boozer_map;
extern bool ullmann_map;

int escape(MAP *m, bool condition(double *x), double ic_box_size)
{	
    FILE *out1 = fopen(filename1, "w");

	int number_of_escape_orbits;
    double x0[2], x1[2];
	double x0_min, x0_max, y0_min, y0_max;
	double deltax, deltay;
	double coordinate, velocity;

	x0_min = m->x0_center - ic_box_size; 
    x0_max = m->x0_center + ic_box_size; 
    y0_min = m->y0_center - ic_box_size; 
    y0_max = m->y0_center + ic_box_size;
	
    deltax = fabs(x0_max - x0_min) / ((double)nx - 1);
    deltay = fabs(y0_max - y0_min) / ((double)ny - 1);

	number_of_escape_orbits = 0;

	coordinate = x0_min; 
    for (int i = 0; i < nx; i++)
    {
    	//printf("Calculating set %d of %d\n", i + 1, nx);
		velocity = y0_min;
        for (int j = 0; j < ny; j++)
        {
			// print progress on terminal
	        //print_prog((double)(j + 1) / (double)ny);

			x0[0] = coordinate;
			x0[1] = velocity;
  
            for (int k = 0; k < iter; k++)
			{	
				if (condition(&x0[1]))
                {
					//printf("test\n");
					number_of_escape_orbits++;
                    fprintf(out1, "%d %1.12e\n", k, x1[0]);
					break;
                }
				// Calling the map .c
				if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
				if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
				if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
				x0[0] = x1[0];
				x0[1] = x1[1];
			}
			// update velocity
			if (ny > 1)
			{
				velocity += deltay;
			}
        }

		// update coordinate
		if (nx > 1)
		{
			coordinate += deltax;
		}
		//printf("\n");
	}

	printf("percentage of escape = %d%%\n", (100 * number_of_escape_orbits)/(nx*ny) );  

	fclose(out1);
    
    return 0;
}