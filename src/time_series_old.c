// --- Function Time Series ---
// Evolve a i.c. and save the time series evolution of a specific coodinate
# include "dynamics.h"

extern char filename1[T], filename2[T];
extern int iter; 
extern bool standard_map;
extern bool boozer_map;
extern bool ullmann_map;
extern bool fermi_ulam_map;

int time_series(MAP *m, int ti, int n_min, int n_max)
{

    FILE *out1 = fopen(filename1, "w");
    FILE *out2 = fopen(filename2, "w");
    
	double x0[2], x1[2];
    double **orbit; 
    double d_x;
    double aux1, aux2;
    double distance_time_series[iter-1];

    alloc_2d_double(&orbit,iter,2);

    x0[0] = m->x0_center;
    x0[1] = m->y0_center;

    fprintf(out2, "%1.12e %1.12e", m->x0_center, m->y0_center);

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
        orbit[k][0] = x1[0];
        orbit[k][1] = x1[1];
	}

    
    
    for (int i = 0; i < iter - 1; i++)
    {	
        d_x = (orbit[i][0] - orbit[i+1][0]);
        if (fabs(d_x) < M_PI) aux1 = d_x;
        else if (d_x > 0) aux1 = 2.0 * M_PI - d_x;
        else aux1 = 2.0 * M_PI + d_x;
        aux1 = aux1 * aux1;
        aux2 = (orbit[i][1] - orbit[i+1][1]);
        aux2 = aux2 * aux2;
        distance_time_series[i] = sqrt(aux1+aux2);
        //printf("%e\n", distance_time_series[i]);
    }

    
    if (n_min != 0 || n_max != 0)
    {
        for (int i = 0; i < iter - 1; i++)
        {
            if (i >= n_min && i <= n_max )
            {
                //if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
		        //&& x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
                //{
                    fprintf(out1, "%1.12e %1.12e %d\n", orbit[i][0], orbit[i][1], i);
                //}
            }
        }
    }
    else
    {
        for (int i = 0; i < iter - 1; i++)
        {
            //if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
		    //&& x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
            //{
                fprintf(out1, "%1.12e %1.12e %d\n", orbit[i][0], orbit[i][1], i);  
            //} 
            //fprintf(out1, "%1.12e %1.12e %d\n", orbit[i][0], orbit[i][1], i);
        }
    }


    printf("\n");

     dealloc_2d_double(&orbit,iter);

	fclose(out1);
    fclose(out2);
    
    return 0;
}