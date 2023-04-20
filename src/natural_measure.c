# include "dynamics.h"

extern char filename1[T];
extern int iter;
extern int nx;
extern int ny;
extern bool standard_map;
extern bool boozer_map;
extern bool ullmann_map;

int natural_measure(MAP *m, int size, double ic_box_size)
{
	FILE *out1 = fopen(filename1, "w");

    int **nm; 
	int inm, jnm;
    double x0[2], x1[2];
	double x0_min, x0_max, y0_min, y0_max;
    double xnm, ynm;
	double coordinate, velocity;
    double dx, dy;

	alloc_2d_int(&nm, size, size);

	x0_min = m->x0_center - ic_box_size; // initial condition minimun x-value
    x0_max = m->x0_center + ic_box_size; // initial condition maximun x-value
    y0_min = m->y0_center - ic_box_size; // initial condition minimun y-value
    y0_max = m->y0_center + ic_box_size; // initial condition maximun y-value
	
    dx = fabs(x0_max - x0_min) / ((double)nx - 1);
    dy = fabs(y0_max - y0_min) / ((double)ny - 1);

    for (int i = 0; i < size; i++)
	{	
		for (int j = 0; j < size; j++)
		{
			nm[i][j] = 0;
		}
	}
    
	coordinate = x0_min; 
    for (int i = 0; i < nx; i++)
    {
		velocity = y0_min; 
		//printf("Calculating set %d of %d\n", i + 1, nx);
        for (int j = 0; j < ny; j++)
        {
			//print_prog((double)(j + 1) / (double)ny);
			x0[0] = coordinate;
			x0[1] = velocity;
	        for (int i = 0; i < iter; i++)
            {
				if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
				if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
				if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
            	x0[0] = x1[0];
				x0[1] = x1[1];
				inm = (int) ((x0[0] - m->x_space_min) * size / (m->x_space_max - m->x_space_min));
                jnm = (int) ((x0[1] - m->y_space_min) * size / (m->y_space_max - m->y_space_min));
				//inm = (int) ((x0[0]-m->x_space_min)/(m->x_space_max-m->x_space_min) * (double)size);
		    	//jnm = (int) ((x0[1]-m->y_space_min)/(m->y_space_max-m->y_space_min) * (double)size);
			
				if (inm >= 0 && inm < size &&
					jnm >= 0 && jnm < size)
                {
						nm[inm][jnm] ++;
				}
           }

            if (ny > 1)
		    {
		    	velocity += dy;
		    }
        }

		if (nx > 1)
		{
		    coordinate += dx;
		}
		//printf("\n");
    }

	for (int i = 0; i < size; i++)
	{
        xnm = (double)i * (m->x_space_max - m->x_space_min) / (double)size + m->x_space_min;
		for (int j = 0; j < size; j++)
		{
            ynm = (double)j * (m->y_space_max - m->y_space_min) / (double)size + m->y_space_min;
			fprintf(out1, "%1.15e %1.15e %d\n", xnm, ynm, nm[i][j]);
		}
		fprintf(out1, "\n");
	}
	
    fclose(out1);

	dealloc_2d_int(&nm, size);
    
    return 0;
}