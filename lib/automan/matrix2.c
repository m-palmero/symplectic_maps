/* -----------------------------------------------------------------
   Allocate memory for vector v(n)
   ----------------------------------------------------------------- */
double*
vector_alloc (int n){
    double *v = malloc (n*sizeof (double));
    return v;
}

/* -----------------------------------------------------------------
   Copy entries of vector v to vector u
   ----------------------------------------------------------------- */
void
vector_copy (double *u, double *v, int dim){
    for (int i=0; i<dim; i++) u[i] = v[i];
}

/* -----------------------------------------------------------------
   The dot product of two vectors
   ----------------------------------------------------------------- */
double
vector_dot (double *u, double *v, int dim){
    double dot = 0.0;
    for (int i=0; i<dim; i++) dot += u[i]*v[i];
    return dot;
}

/* -----------------------------------------------------------------
   The euclidean norm of a vector
   ----------------------------------------------------------------- */
double
vector_norm (double *v, int dim){
    return sqrt (vector_dot (v, v, dim));
}

/* -----------------------------------------------------------------
   Allocate memory for matrix A(m,n)
   ----------------------------------------------------------------- */
double**
matrix_alloc (int m, int n){
    double **A = malloc (m*sizeof (double*));
    for (int i=0; i<m; i++) A[i] = malloc (n*sizeof (double));
    return A;
}

/* -----------------------------------------------------------------
   Copy entries of matrix B(m,n) to A(m,n)
   ----------------------------------------------------------------- */
void
matrix_copy (double **A, double **B, int m, int n){
    for (int i=0; i<m; i++) vector_copy (A[i], B[i], n);
}

/* -----------------------------------------------------------------
   Matrix product AB(m,n) = A(m,l) x B(l,n)
   ----------------------------------------------------------------- */
void
matrix_product (double **AB, double **A, double **B,
		int m, int l, int n){
    int i, j, k;
    for (i=0; i<m; i++){
	for (j=0; j<n; j++){
	    AB[i][j] = 0.0;
	    for (k=0; k<l; k++){
		AB[i][j] += A[i][k]*B[k][j];
	    }
	}
    }
}

/* -----------------------------------------------------------------
   Free memory for matrix A(m,n) : does it require a pointer ***A?
   ----------------------------------------------------------------- */
void matrix_free (double **A, int n){
    for (int i=0; i<n; i++) free (A[i]);
    free (A);
}

/* -------------------------------------------------
   Returns the eigenvalues and eigenvectors of A[2][2]
   in ascending order of the eigenvalue magnitude. The
   eigenvectors are set in the quadrants 1 or 3.
   ------------------------------------------------- */
void
eigensystem (double **A, double **V, double *v){
  int i;
  double a, b, c, d, lp, lm, discriminant, l0, l1;
  
  a = A[0][0]; b = A[0][1]; c = A[1][0]; d = A[1][1];
  discriminant = (a + d)*(a + d) + 4.0*(b*c - a*d);
  if (discriminant > 0){
    lp = 0.5*(a + d + sqrt(discriminant));
    lm = 0.5*(a + d - sqrt(discriminant));
    if (lp*lm > 0.0){
      fprintf (stderr, "Warning! eigenvalues with same sign.\n");
    }
    if (fabs(lp) > fabs(lm)){ v[0] = lm; v[1] = lp; }
    else { v[0] = lp; v[1] = lm; }
    V[0][0] = 1.0;
    V[0][1] = c/(v[0]-d);
    V[1][0] = 1.0;
    V[1][1] = c/(v[1]-d);
    /* Normalize eigenvectors */
    l0 = vector_norm (V[0], 2); l1 = vector_norm (V[1], 2);
    for (i=0; i<2; i++){
      V[0][i] /= l0;
      V[1][i] /= l1;
    }
  }
  else{
    fprintf (stderr, "Complex eigenvalues, leaving now...\n");
  }
}
