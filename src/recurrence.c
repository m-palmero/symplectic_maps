# include "symplectic_maps_analysis.h"

extern char data_file1[T], data_file2[T], data_file3[T], data_file4[T], data_file5[T], data_file6[T];
extern int iter;
extern int nx;
extern int ny;
extern bool standard_map;
extern bool boozer_map;
extern bool ullmann_map;
extern bool fermi_ulam_map;

int recurrence_rate_ensemble(MAP *m, double ensemble_size, int iter_effect_size, double multi, double eps)
{
	FILE *out1 = fopen(data_file1, "w");
	FILE *out2 = fopen(data_file2, "w");
	FILE *out3 = fopen(data_file3, "w");
	FILE *out4 = fopen(data_file4, "w");
	FILE *out5 = fopen(data_file5, "w");

	double **orbit;
	double aux1, aux2, d_x, dist;
	double x0[2], x1[2];
	double x, y;
    double x0_min, x0_max, y0_min, y0_max;
	double deltax, deltay;
	double variance, variance_square;
	double aux_sum_variance, aux_standard_deviation;
	double average_rec_rate;
	double standard_deviation;
	double lower_limit;
	double upper_limit;

//Ensemble of ICs:
    x0_min = m->x0_center - ensemble_size; 
    x0_max = m->x0_center + ensemble_size;
    deltax = fabs(x0_max - x0_min) / ((double)nx - 1);
    y0_min = m->y0_center - ensemble_size; 
    y0_max = m->y0_center + ensemble_size;
    deltay = fabs(y0_max - y0_min) / ((double)ny - 1);

    bool ensemble_in_x, ensemble_in_y;
    if (nx > 1) ensemble_in_x = true;
    else ensemble_in_x = false;
    if (ny > 1) ensemble_in_y = true;
    else ensemble_in_y = false;

    if (ensemble_in_x)
    {
        double rec_rate[nx];
        y = m->y0_center;
        x = x0_min;

        printf("Progress:\n");
        for (int l = 0; l < nx; l++)
        {
            print_prog((double)(l + 1) / (double)nx);
            x0[0] = x;
            x0[1] = y;
            int count = 0;
            for (int k = 0; k < iter; k++)
            {
                if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
                if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
                if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
                if (fermi_ulam_map) sfu_map_eqs(x0, x1, &m->parameter);
                x0[0] = x1[0];
                x0[1] = x1[1];
                if (x1[1] > 1.0) break;
                else
                {
                    if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
                    && x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
                    {
                        if (count == 0)
                        {
                            alloc_2d_double(&orbit, 1, 2);
                        }
                        else
                        {
                            orbit = realloc(orbit, (count + 1) * sizeof(double*));
                            orbit[count] = malloc(2 * sizeof(double));
                        }
                        //NORMALIZATION
                        orbit[count][0] = fabs(x0[0] - m->x_space_min) / (m->x_space_max - m->x_space_min);
                        orbit[count][1] = fabs(x0[1] - m->y_space_min) / (m->y_space_max - m->y_space_min);
                        count++;
                        if (count > iter_effect_size) break;
                    }
                }
            }

            int number_of_recurrences = 0;
            for (int i = 0; i < count; i++)
            {
                for (int j = 0; j < count; j++)
                {
                    aux1 = (orbit[i][0] - orbit[j][0]);
                    aux1 = aux1 * aux1;
                    aux2 = (orbit[i][1] - orbit[j][1]);
                    aux2 = aux2 * aux2;
                    dist = sqrt(aux1 + aux2);
                    if (dist < eps) number_of_recurrences++;
                    rec_rate[l] = (double)number_of_recurrences / (double)(count * count);
                    rec_rate[l] = rec_rate[l] * 100;
                }
            }
            
            //printf("RR = %1.6f\n", rec_rate[l]);

            fprintf(out1, "%1.16e %1.6f \n", x - m->x0_center, rec_rate[l]);

            average_rec_rate += rec_rate[l] / (double)(nx * ny);

            if (nx > 1) x += deltax;

            dealloc_2d_double(&orbit, count);
        }

        printf("\n");
        
        for (int l = 0; l < nx; l++)
        {
            variance = rec_rate[l] - average_rec_rate;
            variance_square = variance * variance;
            aux_sum_variance += variance_square ;
            aux_standard_deviation = aux_sum_variance / (double) (nx - 1);
            standard_deviation = sqrt(aux_standard_deviation);
            lower_limit = average_rec_rate - multi * standard_deviation;
            upper_limit = average_rec_rate + multi * standard_deviation;
        }

        int number_of_peaks = 0;
        for (int l = 0; l < nx; l++)
        {
            //if (rec_rate[l] > upper_limit || rec_rate[l] < lower_limit)
            if (rec_rate[l] > upper_limit)
            {
                number_of_peaks++;
            }  
        }

        printf("%d peaks for %1.4f\n",number_of_peaks, m->parameter);

        fprintf(out2, " AVG_RR = %1.10f\n STD_RR = %1.10f\n Upper limit = %1.10f\n Lower limit = %1.10f\n Number of peaks = %d", average_rec_rate, standard_deviation, upper_limit, lower_limit, number_of_peaks);

        for (double x = x0_min; x <= x0_max; x += deltax)
        {
            fprintf(out3, "%1.16e %1.6f \n", x - m->x0_center, average_rec_rate);
            fprintf(out4, "%1.16e %1.6f \n", x - m->x0_center, lower_limit);
            fprintf(out5, "%1.16e %1.6f \n", x - m->x0_center, upper_limit);
        }

        fclose(out1);
        fclose(out2);
        fclose(out3);
        fclose(out4);
        fclose(out5);
    }

    if (ensemble_in_y)
    {
        double rec_rate[ny];
        y = y0_min;
        x = m->x0_center;

        printf("Progress:\n");
        for (int l = 0; l < ny; l++)
        {
            print_prog((double)(l + 1) / (double)ny);
            x0[0] = x;
            x0[1] = y;
            int count = 0;
            for (int k = 0; k < iter; k++)
            {
                if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
                if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
                if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
                if (fermi_ulam_map) sfu_map_eqs(x0, x1, &m->parameter);
                x0[0] = x1[0];
                x0[1] = x1[1];
                if (boozer_map)
                {
                    if (x1[1] > 1.0) break;
                    else
                    {
                        if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
                        && x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
                        {
                            if (count == 0)
                            {
                                alloc_2d_double(&orbit, 1, 2);
                            }
                            else
                            {
                                orbit = realloc(orbit, (count + 1) * sizeof(double*));
                                orbit[count] = malloc(2 * sizeof(double));
                            }
                            //NORMALIZATION
                            orbit[count][0] = fabs(x0[0] - m->x_space_min) / (m->x_space_max - m->x_space_min);
                            orbit[count][1] = fabs(x0[1] - m->y_space_min) / (m->y_space_max - m->y_space_min);
                            count++;
                            if (count > iter_effect_size) break;
                        }
                    }
                }
                if (ullmann_map)
                {
                    if (x1[1] < 0.0) break;
                    else
                    {
                        if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
                        && x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
                        {
                            if (count == 0)
                            {
                                alloc_2d_double(&orbit, 1, 2);
                            }
                            else
                            {
                                orbit = realloc(orbit, (count + 1) * sizeof(double*));
                                orbit[count] = malloc(2 * sizeof(double));
                            }
                            //NORMALIZATION
                            orbit[count][0] = fabs(x0[0] - m->x_space_min) / (m->x_space_max - m->x_space_min);
                            orbit[count][1] = fabs(x0[1] - m->y_space_min) / (m->y_space_max - m->y_space_min);
                            count++;
                            if (count > iter_effect_size) break;
                        }
                    }
                }
                
            }

            int number_of_recurrences = 0;
            for (int i = 0; i < count; i++)
            {
                for (int j = 0; j < count; j++)
                {
                    aux1 = (orbit[i][0] - orbit[j][0]);
                    aux1 = aux1 * aux1;
                    aux2 = (orbit[i][1] - orbit[j][1]);
                    aux2 = aux2 * aux2;
                    dist = sqrt(aux1 + aux2);
                    if (dist < eps) number_of_recurrences++;
                    rec_rate[l] = (double)number_of_recurrences / (double)(count * count);
                    rec_rate[l] = rec_rate[l] * 100;
                }
            }

            fprintf(out1, "%1.16e %1.6f \n", y - m->y0_center, rec_rate[l]);

            average_rec_rate += rec_rate[l] / (double)(nx * ny);

            if (ny > 1) y += deltay;

            dealloc_2d_double(&orbit, count);
        }

        printf("\n");
        
        for (int l = 0; l < ny; l++)
        {
            variance = rec_rate[l] - average_rec_rate;
            variance_square = variance * variance;
            aux_sum_variance += variance_square ;
            aux_standard_deviation = aux_sum_variance / (double) (ny - 1);
            standard_deviation = sqrt(aux_standard_deviation);
            lower_limit = average_rec_rate - multi * standard_deviation;
            upper_limit = average_rec_rate + multi * standard_deviation;
        }

        int number_of_peaks = 0;
        for (int l = 0; l < ny; l++)
        {
            //if (rec_rate[l] > upper_limit || rec_rate[l] < lower_limit)
            if (rec_rate[l] > upper_limit)
            {
                number_of_peaks++;
            }  
        }

        printf("%d peaks for %1.4f\n",number_of_peaks, m->parameter);

        fprintf(out2, " AVG_RR = %1.10f\n STD_RR = %1.10f\n Upper limit = %1.10f\n Lower limit = %1.10f\n Number of peaks = %d", average_rec_rate, standard_deviation, upper_limit, lower_limit, number_of_peaks);

        for (double y = y0_min; y <= y0_max; y += deltay)
        {
            fprintf(out3, "%1.16e %1.6f \n", y - m->y0_center, average_rec_rate);
            fprintf(out4, "%1.16e %1.6f \n", y - m->y0_center, lower_limit);
            fprintf(out5, "%1.16e %1.6f \n", y - m->y0_center, upper_limit);
        }

        fclose(out1);
        fclose(out2);
        fclose(out3);
        fclose(out4);
        fclose(out5);
    }

    return 0;
}


int find_sticky_orbits(MAP *m, double ensemble_size, int iter_effect_size, double multi, double eps)
{
	FILE *out1 = fopen(data_file1, "w");

	double **orbit;
	double aux1, aux2, d_x, dist;
	double x0[2], x1[2];
	double x, y;
    double x0_min, x0_max, y0_min, y0_max;
	double deltax, deltay;
	double variance, variance_square;
	double aux_sum_variance, aux_standard_deviation;
	double average_rec_rate;
	double standard_deviation;
	double lower_limit;
	double upper_limit;

//Ensemble of ICs:
    x0_min = m->x0_center - ensemble_size; 
    x0_max = m->x0_center + ensemble_size;
    deltax = fabs(x0_max - x0_min) / ((double)nx - 1);
    y0_min = m->y0_center - ensemble_size; 
    y0_max = m->y0_center + ensemble_size;
    deltay = fabs(y0_max - y0_min) / ((double)ny - 1);

    bool ensemble_in_x, ensemble_in_y;
    if (nx > 1) ensemble_in_x = true;
    else ensemble_in_x = false;
    if (ny > 1) ensemble_in_y = true;
    else ensemble_in_y = false;

    if (ensemble_in_x)
    {
        double rec_rate[nx];
        y = m->y0_center;
        x = x0_min;

        printf("Progress:\n");
        for (int l = 0; l < nx; l++)
        {
            print_prog((double)(l + 1) / (double)nx);
            x0[0] = x;
            x0[1] = y;
            int count = 0;
            for (int k = 0; k < iter; k++)
            {
                if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
                if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
                if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
                if (fermi_ulam_map) sfu_map_eqs(x0, x1, &m->parameter);
                x0[0] = x1[0];
                x0[1] = x1[1];
                if (x1[1] > 1.0) break;
                else
                {
                    if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
                    && x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
                    {
                        if (count == 0)
                        {
                            alloc_2d_double(&orbit, 1, 2);
                        }
                        else
                        {
                            orbit = realloc(orbit, (count + 1) * sizeof(double*));
                            orbit[count] = malloc(2 * sizeof(double));
                        }
                        //NORMALIZATION
                        orbit[count][0] = fabs(x0[0] - m->x_space_min) / (m->x_space_max - m->x_space_min);
                        orbit[count][1] = fabs(x0[1] - m->y_space_min) / (m->y_space_max - m->y_space_min);
                        count++;
                        if (count > iter_effect_size) break;
                    }
                }
            }

            int number_of_recurrences = 0;
            for (int i = 0; i < count; i++)
            {
                for (int j = 0; j < count; j++)
                {
                    aux1 = (orbit[i][0] - orbit[j][0]);
                    aux1 = aux1 * aux1;
                    aux2 = (orbit[i][1] - orbit[j][1]);
                    aux2 = aux2 * aux2;
                    dist = sqrt(aux1 + aux2);
                    if (dist < eps) number_of_recurrences++;
                    rec_rate[l] = (double)number_of_recurrences / (double)(count * count);
                    rec_rate[l] = rec_rate[l] * 100;
                }
            }
            
            average_rec_rate += rec_rate[l] / (double)(nx * ny);

            if (nx > 1) x += deltax;

            dealloc_2d_double(&orbit, count);
        }
        
        
        for (int l = 0; l < nx; l++)
        {
            variance = rec_rate[l] - average_rec_rate;
            variance_square = variance * variance;
            aux_sum_variance += variance_square ;
            aux_standard_deviation = aux_sum_variance / (double) (nx - 1);
            standard_deviation = sqrt(aux_standard_deviation);
            lower_limit = average_rec_rate - multi * standard_deviation;
            upper_limit = average_rec_rate + multi * standard_deviation;
        }
        
        int c = 0;
        int *f;
        for (int l = 0; l < nx; l++)
        {
            if (rec_rate[l] > upper_limit)
            {
                //printf("stiky at l = %d; RR = %1.6f\n",l, rec_rate[l]);
                if (c == 0)
                {
                    f = (int*)calloc(1, sizeof(int));
                    //alloc_1d_int(&f, 1);
                }
                else
                {
                    f = realloc(f, (c + 1) * sizeof(int));
                    //f = (int *) malloc(sizeof(int));
                }
                f[c] = l;
                c++;
            }    
        }

        double sx0[2],sx1[2];
        for (int k = 0; k < c; k++)
        {
            y = m->y0_center;
            x = x0_min;
            for (int l = 0; l < nx; l++)
            {
                if(l == f[k])
                {   
                    //printf("stiky at l = %d; RR = %1.6f\n",l, rec_rate[l]);
                    sx0[0] = x;
                    sx0[1] = y;
                    for (int n = 0; n < iter; n++)
                    {
                        if (boozer_map) boozer_map_eqs(sx0, sx1, &m->parameter);
                        if (ullmann_map) ullmann_map_eqs(sx0, sx1, &m->parameter);
                        sx0[0] = sx1[0];
                        sx0[1] = sx1[1];
                        if(sx1[0] <= m->x_space_max && sx1[0] >= m->x_space_min 
                        && sx1[1] <= m->y_space_max && sx1[1] >= m->y_space_min)
                        {
                            fprintf(out1, "%1.16e %1.16e %1.6f %d\n", sx0[0], sx0[1], rec_rate[l], n);
                        }
                    }
                }
                if (nx > 1) x += deltax;
            }
        }

        free(f);

        printf("\n");

        fclose(out1);
    }
    
    if (ensemble_in_y)
    {
        double rec_rate[ny];
        y = y0_min;
        x = m->x0_center;

        printf("Progress:\n");
        for (int l = 0; l < ny; l++)
        {
            print_prog((double)(l + 1) / (double)ny);
            x0[0] = x;
            x0[1] = y;
            int count = 0;
            for (int k = 0; k < iter; k++)
            {
                if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
                if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
                if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
                if (fermi_ulam_map) sfu_map_eqs(x0, x1, &m->parameter);
                x0[0] = x1[0];
                x0[1] = x1[1];
                if (x1[1] > 1.0) break;
                else
                {
                    if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
                    && x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
                    {
                        if (count == 0)
                        {
                            alloc_2d_double(&orbit, 1, 2);
                        }
                        else
                        {
                            orbit = realloc(orbit, (count + 1) * sizeof(double*));
                            orbit[count] = malloc(2 * sizeof(double));
                        }
                        //NORMALIZATION
                        orbit[count][0] = fabs(x0[0] - m->x_space_min) / (m->x_space_max - m->x_space_min);
                        orbit[count][1] = fabs(x0[1] - m->y_space_min) / (m->y_space_max - m->y_space_min);
                        count++;
                        if (count > iter_effect_size) break;
                    }
                }
            }

            int number_of_recurrences = 0;
            for (int i = 0; i < count; i++)
            {
                for (int j = 0; j < count; j++)
                {
                    aux1 = (orbit[i][0] - orbit[j][0]);
                    aux1 = aux1 * aux1;
                    aux2 = (orbit[i][1] - orbit[j][1]);
                    aux2 = aux2 * aux2;
                    dist = sqrt(aux1 + aux2);
                    if (dist < eps) number_of_recurrences++;
                    rec_rate[l] = (double)number_of_recurrences / (double)(count * count);
                    rec_rate[l] = rec_rate[l] * 100;
                }
            }

            average_rec_rate += rec_rate[l] / (double)(nx * ny);

            if (ny > 1) y += deltay;

            dealloc_2d_double(&orbit, count);
        }
        
        for (int l = 0; l < ny; l++)
        {
            variance = rec_rate[l] - average_rec_rate;
            variance_square = variance * variance;
            aux_sum_variance += variance_square ;
            aux_standard_deviation = aux_sum_variance / (double) (ny - 1);
            standard_deviation = sqrt(aux_standard_deviation);
            lower_limit = average_rec_rate - multi * standard_deviation;
            upper_limit = average_rec_rate + multi * standard_deviation;  
        }

        
        int c = 0;
        int *f;
        for (int l = 0; l < ny; l++)
        {
            if (rec_rate[l] > upper_limit)
            {
                //printf("stiky at l = %d; RR = %1.6f\n",l, rec_rate[l]);
                if (c == 0)
                {
                    f = (int*)calloc(1, sizeof(int));
                }
                else
                {
                    f = realloc(f, (c + 1) * sizeof(int));
                }
                f[c] = l;
                c++;
            }    
        }

        double sx0[2],sx1[2];
        for (int k = 0; k < c; k++)
        {
            y = y0_min;
            x = m->x0_center;
            for (int l = 0; l < ny; l++)
            {
                if(l == f[k])
                {   
                    sx0[0] = x;
                    sx0[1] = y;
                    for (int n = 0; n < iter; n++)
                    {
                        if (boozer_map) boozer_map_eqs(sx0, sx1, &m->parameter);
                        if (ullmann_map) ullmann_map_eqs(sx0, sx1, &m->parameter);
                        sx0[0] = sx1[0];
                        sx0[1] = sx1[1];
                        if(sx1[0] <= m->x_space_max && sx1[0] >= m->x_space_min 
                        && sx1[1] <= m->y_space_max && sx1[1] >= m->y_space_min)
                        {
                            fprintf(out1, "%1.16e %1.16e %1.6f %d\n", sx0[0], sx0[1], rec_rate[l], n);
                        }
                    }  
                }
                if (ny > 1) y += deltay;
            }
        }
        free(f);

        printf("\n");

        fclose(out1);

    }

    return 0;
    
}



int recurrence_plot(MAP *m, int iter_effect_size, double eps)
{
    FILE *out1 = fopen(data_file1, "w");
    FILE *out2 = fopen(data_file2, "w");

    int **rec_matrix;
    double **orbit;
    double aux1, aux2, dist;
    double x0[2], x1[2];
    double rec_rate;

    x0[0] = m->x0_center;
    x0[1] = m->y0_center;

    int count = 0;
    for (int k = 0; k < iter; k++)
    {
		if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
		if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
		if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
        if (fermi_ulam_map) sfu_map_eqs(x0, x1, &m->parameter);
		x0[0] = x1[0];
		x0[1] = x1[1];
        if (x1[1] > 1.0) break;
        else
        {
            if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
		    && x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
            {
                if (count == 0)
                {
                    alloc_2d_double(&orbit, 1, 2);
                }
                else
                {
                    orbit = realloc(orbit, (count + 1) * sizeof(double*));
                    orbit[count] = malloc(2 * sizeof(double));
                }
                //NORMALIZATION
			    orbit[count][0] = fabs(x0[0] - m->x_space_min) / (m->x_space_max - m->x_space_min);
        	    orbit[count][1] = fabs(x0[1] - m->y_space_min) / (m->y_space_max - m->y_space_min);
                count++;
                if (count > iter_effect_size) break;
		    }
        }

    }

    alloc_2d_int(&rec_matrix, count, count);

    int rec = 0;
    for (int i = 0; i < count; i++)
    {
        for (int j = 0; j < count; j++)
        {
			aux1 = (orbit[i][0] - orbit[j][0]);
			aux1 = aux1 * aux1;
			aux2 = (orbit[i][1] - orbit[j][1]);
			aux2 = aux2 * aux2;
			dist = sqrt(aux1 + aux2);
            if (dist < eps)
			{
				rec_matrix[i][j] = 1;
				rec++;
			}
			else
			{
				rec_matrix[i][j] = 0;
			}
            //fprintf(out2,"\t%d\t", rec_matrix[i=j][j=i]);
            rec_rate = (double)rec / (double)(count * count);
            rec_rate = rec_rate * 100;
        }
    }

    fprintf(out1, "%.7f\n", rec_rate);

    int x_rp[count], y_rp[count];

    for (int i = 0; i < count; i++)
    {
        for (int j = 0; j < count; j++)
        {
            if (rec_matrix[i][j] == 1)
            {
                x_rp[i] = i;
                y_rp[j] = j;
                //fprintf(out2, "%d %d %d\n", x_rp[i], y_rp[j], i);
                if (i >= j)fprintf(out2, "%d %d %d\n", x_rp[i], y_rp[j], i);
            }
        }
    }

    dealloc_2d_int(&rec_matrix, count);
    dealloc_2d_double(&orbit, count);
    
    fclose(out1);
   
    return 0;
}