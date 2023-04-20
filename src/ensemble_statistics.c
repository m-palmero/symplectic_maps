// --- Function Time Series ---
// Evolve a i.c. and save the time series evolution of a specific coodinate
# include "dynamics.h"

extern char filename1[T], filename2[T];
extern int iter;
extern int nx;
extern int ny;
extern bool standard_map;
extern bool boozer_map;
extern bool ullmann_map;
extern bool fermi_ulam_map;

int evolve_ensemble_statistics(MAP *m)
{	
	FILE *out1 = fopen(filename1, "w");
	FILE *out2 = fopen(filename2, "w");

	double x0[2], x1[2];
    double *avg, *avg2, *rog;
    double sum, sum2;
    double x, y, aux1, aux2, aux3, aux4;
    int lapis = 100;

    alloc_1d_double(&avg,iter);
	alloc_1d_double(&avg2,iter);
	alloc_1d_double(&rog,iter);

    for (int k = 0; k < iter; k++) 
    {
        avg[k] = 0.0; 
        avg2[k] = 0.0;
        rog[k] = 0.0;
    }

    y = m->y0_center;

    //fprintf(out2, "%1.12e %1.12e", m->x0_center, m->y0_center);
  
    for (int i = 0; i < nx; i++)
    {
        //printf("Calculating set %d of %d\n", i + 1, nx);
        x = m->x0_center + (2 * M_PI) * (double) i / (double) nx;
        x0[0] = x;
        x0[1] = y;
        sum = 0.0;
        sum2 = 0.0;
        for (int k = 0; k < iter; k++)
        {	
            //print_prog((double)(i + 1) / (double)nx);
            if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
		    if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
		    if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
            if (fermi_ulam_map) sfu_map_eqs(x0, x1, &m->parameter);
            x0[0] = x1[0];
    	    x0[1] = x1[1];
            sum += x0[1];
            sum2 += (x0[1] * x0[1]);
            // Not trajectory averaging
            //avg[k] += x0[1];
            // Averaging on the trajectory
            avg[k] += (sum / (double) k);
            avg2[k] += (sum2 / (double) k);
            rog[k] = sqrt(fabs(avg2[k]-(avg[k]*avg[k])));
	    }
        //printf("\n");
    }

    //printf("%f\n", rog[100]);

    for (int k = 0; k < iter; k++)
    {
        //aux2 = aux2 / (double) nx;
        //printf("%f\n", avg[k]);
        fprintf(out1, "%1.12e\n", avg[k] / (double) nx);
        fprintf(out2, "%1.12e\n", rog[k] / (double) nx);
    }
     
    dealloc_1d_double(&avg);
	dealloc_1d_double(&avg2);
	dealloc_1d_double(&rog);

	fclose(out1);
	fclose(out2);
    
    return 0;
}