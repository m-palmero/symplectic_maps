# include "dynamics.h"

extern char filename1[T], filename2[T], filename3[T];
extern int iter;
extern int nx;
extern int ny;
extern bool standard_map;
extern bool boozer_map;
extern bool ullmann_map;
extern bool fermi_ulam_map;

int network_transitivity(MAP *m, double eps, int window)
{
    FILE *out1 = fopen(filename1, "w"); //Adjencency Matrix
	FILE *out2 = fopen(filename2, "w");
    FILE *out3 = fopen(filename3, "w"); //Transitivity

    int **rec_matrix;
    int **kron_delta;
    int **adj_matrix;
    double *transitivity;
    int x_rp[iter], y_rp[iter];
    double orbit[iter][2]; 
    double aux1, aux2, d_x, dist;
    double x0[2], x1[2];
    double coordinate, velocity;

    alloc_2d_int(&rec_matrix, iter, iter);
	alloc_2d_int(&kron_delta, iter, iter);
    alloc_2d_int(&adj_matrix, iter, iter);
    alloc_1d_double(&transitivity, iter);

    for (int i = 0; i < iter; i++)
	{	
		for (int j = 0; j < iter; j++)
		{
			rec_matrix[i][j] = 0;
            kron_delta[i][j] = 0;
            adj_matrix[i][j] = 0;
		}
	}

    for (int i = 1; i < iter; i++)
	{	
		for (int j = 1; j < iter; j++)
		{ 
            if (i == j)
            {
                kron_delta[i][j] = 1;
            }
            else 
            {
                kron_delta[i][j] = 0;
            }
		}
    }

    coordinate = m->x0_center;
    velocity = m->y0_center;

	x0[0] = coordinate;
	x0[1] = velocity;
   
    for (int n = 0; n < iter; n++)
	{
        if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
		if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
		if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
        if (fermi_ulam_map) sfu_map_eqs(x0, x1, &m->parameter);
        x0[0] = x1[0];
		x0[1] = x1[1];
        //constrainAngle(&x1[0]);
        //if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
		//&& x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
        //{
            orbit[n][0] = x1[0];
            orbit[n][1] = x1[1];
        //}
    }

    for (int i = 0; i < iter; i++)
    {	
        for (int j = 0; j < iter; j++)
        { 
            d_x = (orbit[i][0] - orbit[j][0]);
            if (fabs(d_x) < M_PI) aux1 = d_x;
            else if (d_x > 0) aux1 = 2.0 * M_PI - d_x;
            else aux1 = 2.0 * M_PI + d_x;
            aux1 = aux1 * aux1;
            aux2 = (orbit[i][1] - orbit[j][1]);
            aux2 = aux2 * aux2;
            dist = sqrt(aux1 + aux2);
            if (dist < eps)
            {
                rec_matrix[i][j] = 1;
            }
            else 
            {
                rec_matrix[i][j] = 0;
            }
            adj_matrix[i][j] = rec_matrix[i][j] - kron_delta[i][j];                    
            //fprintf(out1,"\t%d\t", adj_matrix[i][j]);                                    
        }
        //fprintf(out1, "\n");
    }

    int tnum = 0;
    int tden = 0;

    for (int ii = 0; ii < iter; ii++)
    {
        for (int jj = 0; jj < iter; jj++)
        {
            for (int kk = 0; kk < iter; kk++)
            {
                //if (ii != jj && jj != kk)
                //{
                    tnum = tnum + adj_matrix[ii][jj]*adj_matrix[ii][kk]*adj_matrix[jj][kk];
                    tden = tden + adj_matrix[ii][jj]*adj_matrix[ii][kk]*(1-kron_delta[jj][kk]);
                //}
            }
        }
        transitivity[ii] = (double) tnum / (double) tden;
        //sliding_window(&transitivity,iter,100);
        
    }

    for (int i = 0; i < iter; i++)
    {
        transitivity[i] = moving_window(transitivity, iter, window);
        //if (i >= 100)
        fprintf(out3, "%1.16e\n", transitivity[i]); 
    }

    printf("%f\n",transitivity[iter-1]); 
    
    dealloc_2d_int(&rec_matrix, iter);
	dealloc_2d_int(&kron_delta, iter);
    dealloc_2d_int(&adj_matrix, iter);
    dealloc_1d_double(&transitivity);

    fclose(out1);
	fclose(out2);
    fclose(out3);

    return 0;
}
