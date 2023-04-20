// --- Function Time Series ---
// Evolve a i.c. and save the time series evolution of a specific coodinate
# include "dynamics.h"

extern char data_file1[T];
extern int iter; 
extern bool standard_map;
extern bool boozer_map;
extern bool ullmann_map;
extern bool fermi_ulam_map;

int time_series(MAP *m)
{

    FILE *out1 = fopen(data_file1, "w");
    
	double x0[2], x1[2];
    double **orbit; 
    double d_x;
    double aux1, aux2;
    double distance_time_series[iter-1];

    alloc_2d_double(&orbit,iter,2);

    x0[0] = m->x0_center;
    x0[1] = m->y0_center;

    printf("Iteration progress:\n");
  
    for (int k = 0; k < iter; k++)
    {	
        print_prog((double)(k + 1) / (double)iter);
        // Calling the map .c
        if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
		if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
		if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
        if (fermi_ulam_map) sfu_map_eqs(x0, x1, &m->parameter);
        x0[0] = x1[0];
    	x0[1] = x1[1];
        //constrainAngle(&x1[0]);
        if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
		&& x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
        {
            orbit[k][0] = x1[0];
            orbit[k][1] = x1[1];
            fprintf(out1, "%1.12e %1.12e %d\n", orbit[k][0], orbit[k][1], k);
        }      
	}

    printf("\n");

    dealloc_2d_double(&orbit,iter);

	fclose(out1);
    
    return 0;
}