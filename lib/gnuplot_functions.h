# ifndef _gnuplot_functions_h
# define _gnuplot_functions_h

void plot_gnuplot_general(double parameter, int iter, char file_name[]);
void plot_gnuplot_boozer(double parameter, int iter, char filename[]);
void plot_gnuplot_fum(double parameter, int iter, char file_name[]);
void plot_gnuplot_histogram(char file_name[], double parameter);
void plot_gnuplot_time_series(char title[]);
void plot_gnuplot_RP_matrix(char title[]);
void plot_gnuplot_RP_coordinates(char filename[], double parameter);
void plot_gnuplot_RPs_title(char filename[], int n, double eps);


# endif