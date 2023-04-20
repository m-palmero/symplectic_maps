# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <string.h>

#define GNUPLOT "gnuplot -persist"

struct map 
{
    double xi;
    double yi;
    double par;
    double x_space_min;
    double x_space_max; 
    double y_space_min; 
    double y_space_max; 
};

typedef struct map MAP;

void map1(double *x0, double *x1, double *k);
void map2(double *x0, double *x1, double *k);
int evolve_phase_space(MAP *m);
int compute_recurrence_plot(MAP *m, double eps);
int compute_recurrence_plot_from_input(double eps);
int alloc_1d_double(double **x, int n);
int alloc_2d_double(double ***x, int n, int m);
int alloc_2d_int(int ***x, int n, int m);
int dealloc_1d_double(double **x);
int dealloc_2d_double(double ***x, int n);
int dealloc_2d_int(int ***x, int n);
void plot_gnuplot_RP_coordinates();

