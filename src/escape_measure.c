# include "dynamics.h"

extern char data_file1[T], data_file2[T], data_file3[T], data_file4[T];
extern int iter;
extern int nx;
extern int ny;
extern bool standard_map;
extern bool boozer_map;
extern bool ullmann_map;
extern bool fermi_ulam_map;


int escape_measure(MAP *m, int size, double ic_box_size, bool condition(double *x))
{
    FILE *nme = fopen(data_file1, "w");
	FILE *ic = fopen(data_file2, "w");	
	FILE *hist = fopen(data_file3, "w");
	FILE *hist_mu = fopen(data_file4, "w");
	//FILE *trap = fopen(filename5, "w");
	//FILE *comp = fopen(filename6, "w");

	// declare variables
	int i_map, j_map, orbit_size, orbit_size_inside_region;
	int number_of_escape_orbits;
	int **visitation_control;
	int **grid_occupation;
	double x0[2], x1[2];
	double x0_min, x0_max, y0_min, y0_max;
	double coordinate, velocity;
	double box_c_min, box_c_max, box_nc;
	double box_v_min, box_v_max, box_nv;
    double escape_x_min, escape_x_max, escape_y_min, escape_y_max;
	double complexity;
	// double ensemble_sum_test;
	double **measure_matrix;
	double **maximum_escape_time;
	double **visitation_matrix;

	// x0_min = 0;
	// x0_max = 0;
	// y0_min = 1 - 1e-4;
	// y0_max = 1 - 1e-6;

	// x0_min = 0;
	// x0_max = 0;
	// y0_min = 1 - 1e-4;
	// y0_max = 1 - 1e-6;

	x0_min = m->x0_center - ic_box_size; // initial condition minimun x-value
    x0_max = m->x0_center + ic_box_size; // initial condition maximun x-value
    y0_min = m->y0_center - ic_box_size; // initial condition minimun y-value
    y0_max = m->y0_center + ic_box_size; // initial condition maximun y-value

	// define grid prameters
	box_c_min = m->x_space_min;
	box_c_max = m->x_space_max;
	box_v_min = m->y_space_min;
	box_v_max = m->y_space_max;

	// Zoom in special region
	// box_c_min = -0.0004;
	// box_c_max = 0.0004;
	// box_v_min = 0.99755;
	// box_v_max = 0.99820;
	
    box_nc = size;
	box_nv = size;

    // Escape conditions
    // escape_x_min = 0.0;
    // escape_x_max = 2 * M_PI;
    escape_y_min = 0.0;
    // escape_y_max = ;

	// allocate memory to declared arrays (needs "aux_func.c")
	alloc_2d_int(&grid_occupation, box_nc, box_nv);
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
    	printf("Calculating set %d of %d for parameter %.4f\n", i + 1, nx, m->parameter);
		velocity = y0_min; // velocity value
        for (int j = 0; j < ny; j++)
        {
			// print progress on terminal (needs "aux_func.c")
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

			orbit_size = iter; //define orbit size

			orbit_size_inside_region = 0;

			// first orbit evolution to verify if it escapes
			for (int k = 0; k < iter; k++)
			{
				// escape verification
				// if (x0[0] > escape_x_max && x0[0] < escape_x_min &&
                //    x0[1] > escape_y_max && x0[1] < escape_y_min)
				//bool escape_condition = (x0[1] > 1.0);
				//bool escape_condition = (x0[1] < 0.0);
                if (condition(&x0[1]))
                {
					fprintf(hist, "%d\n", k); //print histogram
					number_of_escape_orbits++; // increase escape number
					orbit_size = k; //redifine orbit size
					break;
				}
				else
				{
					i_map = (int) ((x0[0] - box_c_min) * box_nc / (box_c_max - box_c_min));
                	j_map = (int) ((x0[1] - box_v_min) * box_nv / (box_v_max - box_v_min));
				}
                if (i_map >= 0 && i_map < box_nc &&
					j_map >= 0 && j_map < box_nv)
                {
					orbit_size_inside_region++;
				}

				// evolve
				if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
				if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
				if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
				if (fermi_ulam_map) sfu_map_eqs(x0, x1, &m->parameter);

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

					grid_occupation[i_map][j_map] = 1;

					// escape condition
					if (orbit_size < iter)
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
				if (standard_map) standard_map_eqs(x0, x1, &m->parameter);
				if (boozer_map) boozer_map_eqs(x0, x1, &m->parameter);
				if (ullmann_map) ullmann_map_eqs(x0, x1, &m->parameter);
				if (fermi_ulam_map) sfu_map_eqs(x0, x1, &m->parameter);
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

	//fprintf(trap, "%d\n", nx*ny-number_of_escape_orbits); // write number of trapped orbits 

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
			complexity += (maximum_escape_time[i][j] / (double) iter) * (visitation_matrix[i][j] / (double) (number_of_escape_orbits)) * (measure_matrix[i][j] / (double) (nx*ny));

			// write measure
			fprintf(nme, "%1.15e %1.15e %1.15e\n", coordinate, velocity, measure_matrix[i][j] / (double) (nx*ny));

			
			if(grid_occupation[i][j] == 1)
			{
			// // 	// write histogram of measure OBS: Not including zero valued boxes
				fprintf(hist_mu, "%1.15e\n", measure_matrix[i][j] / (double) (nx*ny));
			}

			velocity += fabs(box_v_max - box_v_min) / (double)(box_nv - 1);
        }

		fprintf(nme, "\n");

		coordinate += fabs(box_c_max - box_c_min) / (double)(box_nc - 1);
    }

	//fprintf(comp, "%1.15e\n", complexity); // write complexity 
	printf("Complexity = %.7f\n", complexity);
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

	// deallocate memory of arrays (needs "aux_func.c")
	dealloc_2d_int(&grid_occupation, box_nc);
	dealloc_2d_int(&visitation_control, box_nc);
	dealloc_2d_double(&measure_matrix, box_nc);
	dealloc_2d_double(&visitation_matrix, box_nc);
	dealloc_2d_double(&maximum_escape_time, box_nc);

	// close exit files
    fclose(nme);
	fclose(ic);
	fclose(hist);
	fclose(hist_mu);
	//fclose(trap);
	//fclose(comp);
	// fclose(gro);
	

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