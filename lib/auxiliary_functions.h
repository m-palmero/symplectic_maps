# ifndef _aux_func_h
# define _aux_func_h


void print_prog(double percentage);
void copy (double *x, double *y, int dim);
int alloc_1d_double(double **x, int n);
int alloc_2d_double(double ***x, int n, int m);
int alloc_2d_int(int ***x, int n, int m);
int alloc_1d_int(int **x, int n);
int dealloc_1d_double(double **x);
int dealloc_1d_int(int **x);
int dealloc_2d_double(double ***x, int n);
int dealloc_2d_int(int ***x, int n);
double constrainAngle(double *x);
double moving_window(double *x, int n, int m);


# endif