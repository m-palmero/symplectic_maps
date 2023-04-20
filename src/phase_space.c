// --- Function Phase Space ---
// Evolve an ensamble of i.c. and save a specific region of the space
# include "dynamics.h"

extern char data_file1[T], data_file2[T];
extern int iter;
extern int nx;
extern int ny;
extern bool standard_map;
extern bool boozer_map;
extern bool ullmann_map;
extern bool fermi_ulam_map;

int phase_space(MAP *m, bool localized_ensemble, bool normalized_phase_space, double ic_box_size, int lapis)
{

	FILE *out1 = fopen(data_file1, "w");
	FILE *out2 = fopen(data_file2, "w");

	double x0[2], x1[2];
	double x0_min, x0_max, y0_min, y0_max;
	double deltax, deltay;
	double coordinate, velocity;

	if (nx > 1 || ny > 1)
    {  
		if (localized_ensemble)
		{
			x0_min = m->x0_center - ic_box_size; 
    		x0_max = m->x0_center + ic_box_size; 
    		y0_min = m->y0_center - ic_box_size; 
    		y0_max = m->y0_center + ic_box_size; 
		}
		else 
		{
			x0_min = m->x0_min;
    		x0_max = m->x0_max;
        	y0_min = m->y0_min;
        	y0_max = m->y0_max;
		}
	}else exit;

    deltax = fabs(x0_max - x0_min) / ((double)nx - 1);
    deltay = fabs(y0_max - y0_min) / ((double)ny - 1);

	// loop over initial conditions 
	coordinate = x0_min; // coordinate value
    for (int i = 0; i < nx; i++)
    {
    	//printf("Calculating set %d of %d\n", i + 1, nx);
		velocity = y0_min; // velocity value
        for (int j = 0; j < ny; j++)
        {
			// print progress on terminal
	        //print_prog((double)(j + 1) / (double)ny);

			x0[0] = coordinate;
			x0[1] = velocity;

			// print initial conditions
			fprintf(out1, "%1.16e %1.16e\n", x0[0], x0[1]);
		
			for (int k = 0; k < iter; k++)
			{	
				// Calling the map .c
				if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
				if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
				if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
				if (fermi_ulam_map) sfu_map_eqs(x0, x1, &m->parameter);
				x0[0] = x1[0];
				x0[1] = x1[1];
				// Set the region of the phase space for ploting
                if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
				&& x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
                {
					if (normalized_phase_space)
					{
						double norm_coodinates[2];
						norm_coodinates[0] = fabs(x0[0] - m->x_space_min) / (m->x_space_max - m->x_space_min);
						norm_coodinates[1] = fabs(x0[1] - m->y_space_min) / (m->y_space_max - m->y_space_min);
						if (k%lapis == 0) fprintf(out2, "%1.16e %1.16e\n", norm_coodinates[0], norm_coodinates[1]);
					}
					else
					{
						//if (k%lapis == 0) fprintf(out2, "%1.16e %1.16e\n", x0[0], x0[1]);
						if (k%lapis == 0) fprintf(out2, "%1.16e %1.16e %d\n", x0[0], x0[1], k);
					}
                }
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
		
		fprintf(out1, "\n");
		
		//printf("\n");
	}

	fclose(out1);
	fclose(out2);
    
    return 0;
}