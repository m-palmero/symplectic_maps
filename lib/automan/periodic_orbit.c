/* ----------------------- DEPENDENCES -----------------------------
   # include "map_jacobian_n.c"
   # include "matrix.c"
   # include "gauss_solver.c"
   ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
   MAIN FUNCTION:
   ----------------------------------------------------------------- */
void periodic_orbit (int n, double *x, void *pars);

/* -----------------------------------------------------------------
   A lenveberg step with damping lamb
   ----------------------------------------------------------------- */
void
minimization_step (double lamb, double *x0, double *x1,
		   double **J, void *pars, int n, double *err){
  double **M, **Jt, *y, *rhs, *xmy, *dx;
  Jt = matrix_alloc (2, 2);
  M = matrix_alloc (2, 2);
  y = vector_alloc (2);
  dx = vector_alloc (2);
  rhs = vector_alloc (2);
  
  /* the functional jacobian */
  J[0][0] -= 1.0; J[1][1] -= 1.0;
  matrix_transpose (Jt, J, 2);
  matrix_product (M, Jt, J, 2, 2, 2);
  M[0][0] *= (1.0+lamb); M[1][1] *= (1.0+lamb);
  /* rhs vector */
  map_n (n, x0, y, pars);
  vector_combine (y, 1.0, -1.0, x0, y, 2);
  matrix_product_vector (rhs, Jt, y, 2, 2);
  gauss_solve (M, dx, rhs, 2);
  /* the new point */
  vector_combine (x1, 1.0, 1.0, x0, dx, 2);
  /* the new error */
  map_n (n, x1, y, pars);
  vector_combine (y, 1.0, -1.0, y, x1, 2);
  *err = vector_norm (y, 2);
  
  free (M); free (Jt); free (y);
  free (rhs); free(dx);
}

/* -----------------------------------------------------------------
   A fixed point algorithm
   ----------------------------------------------------------------- */
void
periodic_orbit (int n, double *x, void *pars){
    int count, max_steps;
    double r, lamb, err0, err1, err2, tol;
    double *tmp, *x0, *x1, *x2, *y, **J;
  
    J = matrix_alloc (2, 2);
    x0 = vector_alloc (2);
    x1 = vector_alloc (2);
    x2 = vector_alloc (2);
    y = vector_alloc (2);

    /* the error of the initial condition */
    vector_copy (x0, x, 2);
    map_n (n, x0, y, pars);
    vector_combine (y, 1.0, -1.0, y, x0, 2);
    err0 = vector_norm (y, 2);
  
    /* the Levenberg-Marquardt parameters */
    r = 2.0;
    lamb = 0.01;
    tol = 1e-15;
    max_steps = 100;
  
    count = 0;
    while (err0 > tol && count < max_steps){
	printf ("step = %d  err = %.5e lamb = %.5e\n",
		    count, err0, lamb);
	jacobian_n (J, n, x0, pars);
	minimization_step (lamb, x0, x1, J, pars, n, &err1);
	minimization_step (r*lamb, x0, x2, J, pars, n, &err2);
	if (err1<err0 && err1<err2){
	    err0 = err1;
	    tmp = x0; x0 = x1; x1 = tmp;
	}
	else if (err2<err0 && err2<err1){
	    err0 = err2;
	    tmp = x0; x0 = x2; x2 = tmp;
	    r *= lamb;
	}
	else if (err2 < err1) lamb *= r;
	else lamb /= r;
	count ++;
    }
    printf ("fixed point found after %d steps\n", count);
    printf ("x* = (%.16e  %.16e)\n", x0[0], x0[1]);
    printf ("|map_n(x*)-x*| = %.6e\n", err0);
    /* if max count was reached */
    if (count == max_steps)
	printf ("Max. count reached during fixed point search\n");
    x[0] = x0[0]; x[1] = x0[1]; // update external var
}


