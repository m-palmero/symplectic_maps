
int iteration_space(int n, int nx, int ny, int size, double x0_min, double x0_max, double y0_min, double y0_max, double x_space_min, double x_space_max, double y_space_min, double y_space_max, double *parameter)
{
	FILE *out1 = fopen("results/iter_space_test.dat", "w");

    int nm[size][size], inm, jnm;
    double x0[2], x1[2];
    double xnm, ynm;
	double coordinate, velocity;
    double dx, dy;
	
    dx = fabs(x0_max - x0_min) / ((double)nx - 1);
    dy = fabs(y0_max - y0_min) / ((double)ny - 1);

    for (int i = 0; i < size; i++)
	{	
		for (int j = 0; j < size; j++)
		{
			nm[i][j] = 0;
		}
	}
    
	coordinate = x0_min; 
    for (int i = 0; i < nx; i++)
    {
		velocity = y0_min; 
        for (int j = 0; j < ny; j++)
        {
			x0[0] = coordinate;
			x0[1] = velocity;
	        for (int i = 0; i < n; i++)
            {
			ullmann_map(x0, x1, parameter); 
            x0[0] = x1[0];
			x0[1] = x1[1];
			inm = (int) ((x0[0]-x_space_min)/(x_space_max-x_space_min) * (double)size);
		    jnm = (int) ((x0[1]-y_space_min)/(y_space_max-y_space_min) * (double)size);
			
			if (inm >= 1 && inm <= size)
			{
				if (jnm >= 1 && jnm <= size)
				{
					nm[inm][jnm] = i;
				}
			}
            
           }

            if (ny > 1)
		    {
		    velocity += dy;
		    }
        }

		if (nx > 1)
		{
		    coordinate += dx;
		}
    }

	for (int i = 0; i < size; i++)
	{
        xnm = (double)i * (x_space_max - x_space_min) / (double)size + x_space_min;
		for (int j = 0; j < size; j++)
		{
            ynm = (double)j * (y_space_max - y_space_min) / (double)size + y_space_min;
			fprintf(out1, "%1.15e %1.15e %d\n", xnm, ynm, nm[i][j]);
		}
		fprintf(out1, "\n");
	}
	
    fclose(out1);
    
    return 0;
}
