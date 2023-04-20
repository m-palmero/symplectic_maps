#define PBSTR "=================================================="
#define PBWIDTH 50

#define T 100 //Maximum characteres for file's names.


void print_prog(double percentage)
{
	int val = (int)(percentage * 100);
	int lpad = (int)(percentage * PBWIDTH);
	int rpad = PBWIDTH - lpad;
	printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
	fflush(stdout);
}

void copy (double *x, double *y, int dim)
{
  for (int i=0; i<dim; i++)	x[i] = y[i];
}

int alloc_1d_double(double **x, int n)
{
	*x = (double *)malloc(n * sizeof(double));
	return 0;
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

int alloc_2d_int(int ***x, int n, int m)
{
	*x = (int **)malloc(n * sizeof(int *));
	for (int i = 0; i < n; i++)
	{
		(*x)[i] = (int *)malloc(m * sizeof(int));
	}
	return 0;
}

int alloc_1d_int(int **x, int n)
{
	*x = (int *)malloc(n * sizeof(int));
	return 0;
}


int dealloc_1d_double(double **x)
{
	free(*x);
	return 0;
}

int dealloc_1d_int(int **x)
{
	free(*x);
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

int dealloc_2d_int(int ***x, int n)
{
	for (int i = 0; i < n; i++)
	{
		free((*x)[i]);
	}
	free(*x);
	return 0;
}

double constrainAngle(double *x)
{
    *x = fmod(*x, 2 * M_PI);
    if (*x < 0)
        *x += 2 * M_PI;
    return *x;
}

double moving_window(double *x, int n, int m)
{
	double aux;

	for (int i = 0; i < n - m; i++)
	{
    	aux = 0;
    	for (int j = 0; j < m; j++)
		{
        	aux += x[i+j];
    	}
    	x[i] = aux * (1.0 / m);
	}

	return *x;
}

