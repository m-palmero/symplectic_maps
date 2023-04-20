#include "test.h"

char filename1[100], filename2[100], filename3[100];
bool map1_on = true;
bool map2_on = false; 

MAP m1;
MAP m2;

int n = 1000;

int main()
{
    bool phase_space = true;
    bool recurrence_plot = true;
    bool recurrence_plot_from_input = false;

    if (map1_on)
    {
        m1.xi = -8.2627236730887474e-06;
        m1.yi = -3.1415982862976222e+00;
        m1.par = 1.46;
        m1.x_space_min = 0.0;
        m1.x_space_max = 2 * M_PI;
        m1.y_space_min = 0.0;
        m1.y_space_max = 2 * M_PI;
    } 
    if (map2_on)
    {
        m2.xi = 1.880006567089865e-16;
        m2.yi = 9.971293753820900e-01; 
        m2.par = 0.6;
        m2.x_space_min = -0.005;
        m2.x_space_max = 0.005;
        m2.y_space_min = 0.996;
        m2.y_space_max = 1.0;
    } 

    if(phase_space)
    {
        sprintf(filename1, "ic.dat");
        sprintf(filename2, "ps.dat");
        sprintf(filename3, "ts.dat");

        if (map1_on) evolve_phase_space(&m1);
        if (map2_on) evolve_phase_space(&m2);
    }

    if(recurrence_plot)
    {
        sprintf(filename1, "rp.dat");
        sprintf(filename2, "test.dat");
        double eps = 0.5;

        if (map1_on) compute_recurrence_plot(&m1, eps);
        if (map2_on) compute_recurrence_plot(&m2, eps);

        //plot_gnuplot_RP_coordinates();
    }

}

void map1(double *x0, double *x1, double *k)
{ 
    x1[1] = x0[1] +  *k * sin (x0[0]);
    //x1[1] = fmod(x1[1],2*M_PI);	
    // if (x1[1] > 2*M_PI)
    // {
    //     x1[1] -= 2*M_PI;
    // }
    //     else if (x1[1] < 0)
    // {
    //     x1[1] += 2*M_PI;
    // }

    x1[0] = x0[0] + x1[1] + M_PI;
    //x1[0] = fmod(x1[0],2*M_PI);
//     if (x1[0] > 2*M_PI)
//     {
//         x1[0] -= 2*M_PI;
//     }
//         else if (x1[0] < 0)
//     {
//         x1[0] += 2*M_PI;
//     }
}

void map2(double *x0, double *x1, double *k)
{
    x1[0] = x0[0] - *k * x0[1] * (1 - x0[1]);
	x1[1] = x0[1] + *k * x1[0];
}


int evolve_phase_space(MAP *m)
{
    double x0[2], x1[2];
    double orbit[n][2]; 
    double d_x;
    double aux1, aux2;
    double distance_time_series[n-1];

    FILE *ic = fopen(filename1, "w");
	FILE *ps = fopen(filename2, "w");
    FILE *ts = fopen(filename3, "w");

    x0[0] = m->xi;
    x0[1] = m->yi;
    fprintf(ic, "%1.15e %1.15e\n", x0[0], x0[1]);
    for (int k = 0; k < n; k++)
    {
        if (map1_on) map1(x0,x1,&m->par);
        if (map2_on) map2(x0,x1,&m->par);
	    x0[0] = x1[0];
		x0[1] = x1[1];
        orbit[k][0] = x1[0];
        orbit[k][1] = x1[1];
        if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
		&& x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
        fprintf(ps, "%1.12e %1.12e %d\n", x0[0], x0[1], k);
    }

    for (int i = 0; i < n - 1; i++)
    {	
        d_x = (orbit[i][0] - orbit[i+1][0]);
        if (fabs(d_x) < M_PI) aux1 = d_x;
        else if (d_x > 0) aux1 = 2.0 * M_PI - d_x;
        else aux1 = 2.0 * M_PI + d_x;
        aux1 = aux1 * aux1;
        aux2 = (orbit[i][1] - orbit[i+1][1]);
        aux2 = aux2 * aux2;
        distance_time_series[i] = sqrt(aux1+aux2);
        fprintf(ts, "%1.12e\n", distance_time_series[i]);
    }
    fclose(ic);
    fclose(ps);
    fclose(ts);
}

int compute_recurrence_plot(MAP *m, double eps)
{
    FILE *rp = fopen(filename1, "w");
	//FILE * = fopen(filename2, "w");

    int **rec_matrix;
    int **kron_delta;
    int **adj_matrix;
    int time[n];
    int diag[n];
    int x_rp[n], y_rp[n];
    int lmin;
    int cont, cont1, cont2, tc;
    int sum, prod, dnum, dden;
    double orbit[n][2]; 
    double aux1, aux2, d_x;
    double x0[2], x1[2];
    double coordinate, velocity;
    double aux, dist, rec_rate, det;

    alloc_2d_int(&rec_matrix, n, n);
	alloc_2d_int(&kron_delta, n, n);
    alloc_2d_int(&adj_matrix, n, n);

    for (int i = 0; i < n; i++)
	{	
		for (int j = 0; j < n; j++)
		{
			rec_matrix[i][j] = 0;
            kron_delta[i][j] = 0;
            adj_matrix[i][j] = 0;
		}
	}

    for (int i = 1; i < n; i++)
	{	
		    for (int j = 1; j < n; j++)
		    { 
                if (i == j)
                {
                    kron_delta[i][j] = 1;
                }
                else kron_delta[i][j] = 0;
		    }
    }

	x0[0] =  m->xi;
	x0[1] =  m->yi;
   
    for (int k = 0; k < n; k++)
	{
        if (map1_on) map1(x0, x1, &m->par);
		if (map2_on) map2(x0, x1, &m->par);
        x0[0] = x1[0];
		x0[1] = x1[1];
        //if(x1[0] <= m->x_space_max && x1[0] >= m->x_space_min 
		//&& x1[1] <= m->y_space_max && x1[1] >= m->y_space_min)
        //{
            orbit[k][0] = x1[0];
            orbit[k][1] = x1[1];
        //}
    }

    //Here I start calculating the critical iteration number n_c, but need more work:
    //(Maybe do in a different function)
    cont1 = 0;
    for (int k = 0; k < n; k++)
    {
        //printf("%f\n",fabs(orbit[n+7][0] - orbit[n][0]));
        if (fabs(orbit[k+7][0] - orbit[k][0]) < eps)
        {
            time[k] = k+7;
            cont1 ++; 
            //printf("%d\n",time[n]);   
        }
    }

    cont2 = 0;
    for (int l = 0; l <= cont1; l++)
    {
        //printf("%d\n",time[l+1]);
        if(time[l+1]-time[l]==1)
        {
            cont2++;
        }
    }

    for (int l = 0; l <= cont2; l++)
    {
        //printf("%d\n",time[l]);
        tc = time[cont2];
    }

    //printf("tc = %d\n",tc);  
    //Here ends the n_c calculation    

    cont=0;
    for (int i = 0; i < n; i++)
    {	
        for (int j = 0; j < n; j++)
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
                x_rp[i]=i;
                y_rp[j]=j;
                fprintf(rp,"%d %d\n", x_rp[i], y_rp[j]);
            }
            else rec_matrix[i][j] = 0;
            adj_matrix[i][j] = rec_matrix[i][j] - kron_delta[i][j];      
            rec_rate = (double) cont / (double) (n*n) ;
        }  
    }
    printf("recurrence rate = %1.2f%%\n",rec_rate*100);  

    dealloc_2d_int(&rec_matrix, n);
	dealloc_2d_int(&kron_delta, n);
    dealloc_2d_int(&adj_matrix, n);

    fclose(rp);
	//fclose();

    return 0;
}

void plot_gnuplot_RP_coordinates()
{	
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	
	fprintf(gp, "set terminal pngcairo size 2000,2000 font 'Times,30'\n");
	fprintf(gp, "set size square \n");
	fprintf(gp, "unset key \n");
	//fprintf(gp, "unset xtics \n");
	//fprintf(gp, "unset ytics \n");
	fprintf(gp, "unset colorbox \n");
	fprintf(gp, "set autoscale yfix \n");
	fprintf(gp, "set autoscale xfix \n");
	//fprintf(gp, "set xtics 0,500,1000 \n");
	//fprintf(gp, "set ytics 0,500,1000 \n");
	fprintf(gp, "set palette maxcolors 2\n");
	fprintf(gp, "set palette defined (0 'white', 1 'black')\n");
	fprintf(gp, "set output 'rp.png'\n");
	fprintf(gp, "plot 'rp.dat' w p pt 5 ps 0.8 lc 'black' notitle\n");
	//fprintf(gp, "pause -1");
	fclose(gp);
}

int alloc_2d_double(double ***x, int n, int m)
{
	*x = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++)
	{
		(*x)[i] = (double *)malloc(m * sizeof(double));
	}
	return 0;
}

int dealloc_2d_double(double ***x, int n)
{
	for (int i = 0; i < n; i++)
	{
		free((*x)[i]);
	}
	free(*x);
	return 0;
}

int alloc_1d_double(double **x, int n)
{
	*x = (double *)malloc(n * sizeof(double));
	return 0;
}

int dealloc_1d_double(double **x)
{
	free(*x);
	return 0;
}

int alloc_2d_int(int ***x, int n, int m)
{
	*x = (int **)malloc(n * sizeof(int *));
	for (int i = 0; i < n; i++)
	{
		(*x)[i] = (int *)malloc(m * sizeof(int));
	}
	return 0;
}

int dealloc_2d_int(int ***x, int n)
{
	for (int i = 0; i < n; i++)
	{
		free((*x)[i]);
	}
	free(*x);
	return 0;
}
