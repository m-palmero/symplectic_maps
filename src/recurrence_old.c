# include "dynamics.h"

extern char filename1[T], filename2[T], filename3[T];
extern int iter;
extern int nx;
extern int ny;
extern bool standard_map;
extern bool boozer_map;
extern bool ullmann_map;
extern bool fermi_ulam_map;

int recurrence_plot(MAP *m, double eps_min, double eps_max, double multi, double eps)
{
    FILE *out1 = fopen(filename1, "w");
	FILE *out2 = fopen(filename2, "w");
    FILE *out3 = fopen(filename3, "w");

    int **rec_matrix;
    int time[iter];
    int diag[iter];
    int x_rp[iter], y_rp[iter];
    int lmin;
    int cont, cont1, cont2, tc;
    int sum, prod, dnum, dden;
    double orbit[iter][2]; 
    double aux1, aux2, d_x;
    double x0[2], x1[2];
    double coordinate, velocity;
    double aux, dist, rec_rate, det;

    alloc_2d_int(&rec_matrix, iter, iter);

    for (int i = 0; i < iter; i++)
	{	
		for (int j = 0; j < iter; j++)
		{
			rec_matrix[i][j] = 0;
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

//Here I start calculating the critical iteration number n_c, but needs more work:
//(Maybe do in a different function)
//     cont1 = 0;
//     for (int n = 0; n < iter; n++)
//     {
//         //printf("%f\n",fabs(orbit[n+7][0] - orbit[n][0]));
//         if (fabs(orbit[n+7][0] - orbit[n][0]) < eps)
//         {
//             time[n] = n+7;
//             cont1 ++; 
//             //printf("%d\n",time[n]);   
//         }
//     }

//     cont2 = 0;
//     for (int l = 0; l <= cont1; l++)
//     {
//         //printf("%d\n",time[l+1]);
//         if(time[l+1]-time[l]==1)
//         {
//             cont2++;
//         }
//     }

//     for (int l = 0; l <= cont2; l++)
//     {
//         //printf("%d\n",time[l]);
//         tc = time[cont2];
//     }

//     //printf("tc = %d\n",tc);  
// //Here ends the n_c calculation    

    if (eps_min != 0 || eps_max != 0)
    {
        for (double eps = eps_min ; eps < eps_max ; eps = eps * multi)
        {   
            //printf("%1.8f\n",eps);
            cont=0;
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
                        cont++;
                    }
                    else 
                    {
                        rec_matrix[i][j] = 0;
                    }
                    
                    fprintf(out1,"\t%d\t", rec_matrix[i][j]);
                    
                    rec_rate = (double) cont / (double) (iter*iter) ;
                }
                fprintf(out1, "\n");
            }
            //printf("recurrence rate = %1.2f%%\n",rec_rate*100);
            fprintf(out2,"%1.10f %d %1.10f\n",eps, tc, rec_rate*100);  
        }

    }
    else
    {
        cont=0;
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
                    cont++;
                }
                else 
                {
                    rec_matrix[i][j] = 0;
                }
                    
                fprintf(out1,"\t%d\t", rec_matrix[i][j]);
                    
                rec_rate = (double) cont / (double) (iter*iter) ;
            }
            fprintf(out1, "\n"); 
        }
        printf("recurrence rate = %1.14e\n",rec_rate*100);
    }

    for (int i = 0; i < iter; i++)
    {
        for (int j = 0; j < iter; j++)
        {
            if (rec_matrix[i][j] == 1)
            {
                x_rp[i]=i;
                y_rp[j]=j;
                fprintf(out3,"%d %d\n", x_rp[i], y_rp[j]);
                //printf("%d\n",x_rp[i]);
            }
        }
    }


    dealloc_2d_int(&rec_matrix, iter);

    fclose(out1);
	fclose(out2);
    fclose(out3);

    return 0;
}
