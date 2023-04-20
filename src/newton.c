/* -------------------------------------------------------------
   Dependences /Programs/Matrix/matrix.c
   /Programs/Solvers/gauss_solver.c
   ------------------------------------------------------------- */

/* -------------------------------------------------------------
   The n-dimensional function and its jacobian
   ------------------------------------------------------------- */
void (*newt_f) (double *f, double *x, double *a);
void (*newt_J) (double **J, double *x, double *a);
int NEWT_DIM;        /* dimension of the problem */
int NEWT_MAX_STEPS;  /* maximum search steps */
double NEWT_TOL;     /* tolerance of the root */
double NEWT_REF;     /* refinement parameter */

/* -------------------------------------------------------------
   Setup to assign user-defined functions to internal functions
   ------------------------------------------------------------- */
void newton_setup (void (*f)(double*, double*, double*),
		   void (*J)(double**, double*, double*), int dim){
    /* redefine local functions */
    newt_f = f;
    newt_J = J;
    NEWT_DIM = dim;
    NEWT_TOL = 1e-12; // defaults may be changed after setup
    NEWT_REF = 0.5;
    NEWT_MAX_STEPS = 100;
}

/* -------------------------------------------------------------
   Search the root of newt_f starting at x with params a
   ------------------------------------------------------------- */
void
newton_zero (double *x, double *a){
    int i;
    double tol, eps, r, err0, err1;
    double *tmp, *x0, *x1, *f0, *f1, **J;
    x0 = vector_alloc (NEWT_DIM);
    x1 = vector_alloc (NEWT_DIM);
    f0 = vector_alloc (NEWT_DIM);
    f1 = vector_alloc (NEWT_DIM);
    J = matrix_alloc (NEWT_DIM, NEWT_DIM);
    
    vector_copy (x0, x, NEWT_DIM);
    /* initialize error */
    newt_f (f0, x0, a);
    err0 = vector_dot (f0, f0, NEWT_DIM);
    printf ("err0 = %.5e\n", err0);
    /* damping parameter */
    eps = -1.0;
    for (i=0; i<NEWT_MAX_STEPS && err0 > NEWT_TOL*NEWT_TOL; i++){
	newt_J (J, x0, a);
	gauss_solve (J, x1, f0, NEWT_DIM); // x1 = -dx
	vector_combine (x1, 1.0, eps, x0, x1, NEWT_DIM);
	newt_f (f1, x1, a);
	err1 = vector_dot (f1, f1, NEWT_DIM);
	if (err1 < err0){ // success
	    err0 = err1;
	    tmp = x0; x0 = x1; x1 = tmp;
	    tmp = f0; f0 = f1; f1 = tmp;
	}
	else {
	    eps *= NEWT_REF;
	}
    }
    if (i == NEWT_MAX_STEPS){
	printf ("Max. steps reached during Newton search\n");
    }
    vector_copy (x, x0, NEWT_DIM);
}
