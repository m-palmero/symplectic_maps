#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "aux_func.c"

int calculate_measure_with_exit_condition(int n, int nx, int ny, int size, double x0_min, double x0_max, double y0_min, double y0_max, double *parameter)
{
    char filename[100];
	
	sprintf(filename, "results/natural_measure_exit_%1.3f.dat", *parameter);
	FILE *nme = fopen(filename, "w");	
	sprintf(filename, "results/initial_conditions_%1.3f.dat", *parameter);
	FILE *ic = fopen(filename, "w");	
    sprintf(filename, "results/histogram_%1.3f.dat", *parameter);
	FILE *hist = fopen(filename, "w");
    // sprintf(filename, "results/histogram_measure_%1.3f.dat", *parameter);
	// FILE *hist_mu = fopen(filename, "w");
    sprintf(filename, "results/trapped_orbits_%1.3f.dat", *parameter);
	FILE *trap = fopen(filename, "w");
    sprintf(filename, "results/complexity_%1.3f.dat", *parameter);
	FILE *comp = fopen(filename, "w");
    // sprintf(filename, "results/grid_occupation_%1.3f.dat", *parameter);
	// FILE *gro = fopen(filename, "w");

	// declare variables
	int i_map, j_map, orbit_size, orbit_size_inside_region;
	int number_of_escape_orbits;
	int **visitation_control;
	int **grid_occupation;
	double x0[2], x1[2];
	double coordinate, velocity;
	double box_c_min, box_c_max, box_nc;
	double box_v_min, box_v_max, box_nv;
    double escape_x_min, escape_x_max, escape_y_min, escape_y_max;
	double complexity;
	// double ensemble_sum_test;
	double **measure_matrix;
	double **maximum_escape_time;
	double **visitation_matrix;

	// define grid prameters
	box_c_min = 0.0;
	box_c_max = 2 * M_PI;
	box_v_min = 0;
	box_v_max = 0.4;

	// Zoom in special region
	// box_c_min = -0.0004;
	// box_c_max = 0.0004;
	// box_v_min = 0.99755;
	// box_v_max = 0.99820;
	
    box_nc = size;
	box_nv = size;

    // Escape conditions
    // escape_x_min = 0.0;
    // escape_x_max = DOIS_PI;
    escape_y_min = 0;
    // escape_y_max = ;

	// allocate memory to declared arrays (needs auxiliar_functions_vmo.c)
	// alloc_2d_int(&grid_occupation, box_nc, box_nv);
	alloc_2d_int(&visitation_control, box_nc, box_nv);
	alloc_2d_double(&measure_matrix, box_nc, box_nv);
	alloc_2d_double(&visitation_matrix, box_nc, box_nv);
	alloc_2d_double(&maximum_escape_time, box_nc, box_nv);

	// initialize arrrays
    for (int i = 0; i < box_nc; i++)
    {
        for (int j = 0; j < box_nv; j++)
        {
			// grid_occupation[i][j] = 0;
            measure_matrix[i][j] = 0.0;
			visitation_matrix[i][j] = 0.0;
			maximum_escape_time[i][j] = 0.0;
        }
    }

	number_of_escape_orbits = 0; // initialize calculation of escape orbits number

	// loop over initial conditions 
	coordinate = x0_min; // coordinate value
    for (int i = 0; i < nx; i++)
    {
    	printf("Calculating set %d of %d\n", i + 1, nx);
		velocity = y0_min; // velocity value
        for (int j = 0; j < ny; j++)
        {
			// print progress on terminal (needs auxiliar_functions_vmo.c)
	        print_prog((double)(j + 1) / (double)ny);

			x0[0] = coordinate;
			x0[1] = velocity;

			// print initial conditions
			fprintf(ic, "%1.15e %1.15e\n", x0[0], x0[1]);

			// zero visitation control matrix
			for (int i = 0; i < box_nc; i++)
			{
				for (int j = 0; j < box_nv; j++)
				{
					visitation_control[i][j] = 0;
				}
			}

			orbit_size = n; //define orbit size

			orbit_size_inside_region = 0;

			// first orbit evolution to verify if it escapes
			for (int k = 0; k < n; k++)
			{
				// escape verification
				// if (x0[0] > escape_x_max && x0[0] < escape_x_min &&
                //    x0[1] > escape_y_max && x0[1] < escape_y_min)
                if (x0[1] < escape_y_min)
				{
					fprintf(hist, "%d\n", k); //print histogram
					number_of_escape_orbits++; // increase escape number
					orbit_size = k; //redifine orbit size
					break;
				}

				i_map = (int) ((x0[0] - box_c_min) * box_nc / (box_c_max - box_c_min));
                j_map = (int) ((x0[1] - box_v_min) * box_nv / (box_v_max - box_v_min));

                if (i_map >= 0 && i_map < box_nc &&
					j_map >= 0 && j_map < box_nv)
                {
					orbit_size_inside_region++;
				}

				// evolve
				ullmann_map (x0, x1, parameter);

				// update
    			x0[0] = x1[0];
    			x0[1] = x1[1];
			}

			// restart initial condition
			x0[0] = coordinate;
			x0[1] = velocity;

			// second orbit evolution
            for (int k = 0; k < orbit_size; k++)
            {
				// gridfy orbit using x0
                i_map = (int) ((x0[0] - box_c_min) * box_nc / (box_c_max - box_c_min));
                j_map = (int) ((x0[1] - box_v_min) * box_nv / (box_v_max - box_v_min));

                if (i_map >= 0 && i_map < box_nc &&
					j_map >= 0 && j_map < box_nv)
                {
					// calculate measure
					measure_matrix[i_map][j_map] += 1.0 / (double) orbit_size_inside_region; 

					// grid_occupation[i_map][j_map] = 1;

					// escape condition
					if (orbit_size < n)
					{
						// calculate visitation
						if(visitation_control[i_map][j_map] == 0)
						{
							visitation_matrix[i_map][j_map] += 1.0;
							visitation_control[i_map][j_map] = 1;
						}

						// calculate maximum escape time
						if (orbit_size > maximum_escape_time[i_map][j_map])
						{
							maximum_escape_time[i_map][j_map] = orbit_size;
						}
					}
                }

				// evolve and update
				ullmann_map (x0, x1, parameter);
    			x0[0] = x1[0];
    			x0[1] = x1[1];

            }

			// update velocity
			if (ny > 1)
			{
				velocity += fabs(y0_max - y0_min) / (double)(ny - 1);
			}
        }

		// update coordinate
		if (nx > 1)
		{
			coordinate += fabs(x0_max - x0_min) / (double)(nx - 1);
		}
		
		fprintf(ic, "\n");

		printf("\n");
    }

	fprintf(trap, "%d\n", nx*ny-number_of_escape_orbits); // write number of trapped orbits 

	// define complexity
	complexity = 0.0;

	// loop over grid
	coordinate = box_c_min;
    for (int i = 0; i < box_nc; i++)
    {
		velocity = box_v_min;

        for (int j = 0; j < box_nv; j++)
        {
			// calculate complexity
			complexity += (maximum_escape_time[i][j] / (double) n) * (visitation_matrix[i][j] / (double) (number_of_escape_orbits)) * (measure_matrix[i][j] / (double) (nx*ny));

			// write measure
			fprintf(nme, "%1.15e %1.15e %1.15e\n", coordinate, velocity, measure_matrix[i][j] / (double) (nx*ny));

			
			// if(grid_occupation[i][j] == 1)
			// {
			// 	// write occuppied cells
			// 	fprintf(gro, "%1.15e %1.15e\n", coordinate, velocity);

			// 	// write histogram of measure OBS: Not including zero valued boxes
			// 	fprintf(hist_mu, "%1.15e\n", measure_matrix[i][j] / (double) (nx*ny));
			// }

			velocity += fabs(box_v_max - box_v_min) / (double)(box_nv - 1);
        }

		fprintf(nme, "\n");

		coordinate += fabs(box_c_max - box_c_min) / (double)(box_nc - 1);
    }

	fprintf(comp, "%1.15e\n", complexity); // write complexity 
	
	// ensemble_sum_test = 0.0;

	// // calculate I
	// I[0] = 0.0;
	// I[1] = 0.0;
	// I[2] = 0.0;
    // for (int i = 0; i < box_nc; i++)
    // {
    //     for (int j = 0; j < box_nv; j++)
    //     {
	// 		if(grid_occupation[i][j] == 1)
	// 		{
	// 			I[0] += 1.0;
	// 			I[1] += (measure_matrix[i][j] / (double) (nc*nv)) * log(measure_matrix[i][j] / (double) (nx*ny));
	// 			I[2] += (measure_matrix[i][j] / (double) (nc*nv)) * (measure_matrix[i][j] / (double) (nx*ny));
	// 			ensemble_sum_test += (measure_matrix[i][j] / (double) (nx*ny));
	// 		}
	// 	}
    // }

	// printf("ensemble_sum_test = %1.3e\n", ensemble_sum_test);

	// deallocate memory of arrays (needs auxiliar_functions_vmo.c)
	// dealloc_2d_int(&grid_occupation, box_nc);
	dealloc_2d_int(&visitation_control, box_nc);
	dealloc_2d_double(&measure_matrix, box_nc);
	dealloc_2d_double(&visitation_matrix, box_nc);
	dealloc_2d_double(&maximum_escape_time, box_nc);

	// close exit files
    fclose(nme);
	fclose(ic);
	fclose(hist);
	fclose(trap);
	fclose(comp);
	// fclose(gro);
	// fclose(hist_mu);

    return 0;
}

// int calculate_escape_dimension(void *params, int n, double coordinate_min, double coordinate_max, double velocity_min, double velocity_max, int nc, int nv)
// {
// 	FILE *D0 = fopen("results/box_counting_dimension.dat", "w");
// 	FILE *D1 = fopen("results/information_dimension.dat", "w");
// 	FILE *D2 = fopen("results/correlation_dimension.dat", "w");

// 	int size;
// 	double I[3];

// 	// NO PARALLEL

// 	// size = 128; // 2^7
// 	// for (int i = 0; i < 5; i++)
// 	// {
// 	// 	printf("Calculating set %d of %d\n", i + 1, 5);

// 	// 	calculate_natural_measure_with_exit_condition(params, n, coordinate_min, coordinate_max, velocity_min, velocity_max, nc, nv, size, I);

// 	// 	fprintf(D0, "%1.15e %1.15e\n", log((double) size), log(I[0]));
// 	// 	fprintf(D1, "%1.15e %1.15e\n", log(1.0 / (double) size), I[1]);
// 	// 	fprintf(D2, "%1.15e %1.15e\n", log(1.0 / (double) size), log(I[2]));
	
// 	// 	size *= 2;
// 	// }

// 	// PARALLEL VERSION

// 	#pragma omp parallel private(size, I)
//   	{
// 	#pragma omp for
// 		for (int i = 7; i < 12; i++)
// 		{
// 			printf("Calculating set %d of %d\n", i - 6, 5);

// 			size = (int) pow(2.0, (double) i);

// 			calculate_natural_measure_with_exit_condition(params, n, coordinate_min, coordinate_max, velocity_min, velocity_max, nc, nv, size, I);

// 			#pragma omp critical
// 			{
// 				fprintf(D0, "%1.15e %1.15e\n", log((double) size), log(I[0]));
// 				fprintf(D1, "%1.15e %1.15e\n", log(1.0 / (double) size), I[1]);
// 				fprintf(D2, "%1.15e %1.15e\n", log(1.0 / (double) size), log(I[2]));
// 			}
// 		}
// 	}

// 	fclose(D0);
// 	fclose(D1);
// 	fclose(D2);

// 	return 0;
// }