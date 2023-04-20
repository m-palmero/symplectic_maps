# include "optimization_gco.h"
/* --------------------------------------------------------------------------
    User defined function
   -------------------------------------------------------------------------- */
/* Evalates the function at all locations */
void lev_marq_function (double *f, double *x, double *b, int nx, int nb);
/* The jacobian vectors at every location */
void lev_marq_jacobian (double **J, double *x, double *b, int nx, int nb);

/* --------------------------------------------------------------------------
    Contents
   -------------------------------------------------------------------------- */
   optimizer create_optimizer (int nx, int nb, double *y, double *x,
     double *b);
   void lev_marq_step (optimizer OP, double *err1, double *err2);
   void lev_marq_optimize (optimizer OP, double *err, double derr_min,
     int max_steps);

/* --------------------------------------------------------------------------
    Create fundamental structure for optimization functions
   -------------------------------------------------------------------------- */
optimizer
create_optimizer (int nx, int nb, double *y, double *x, double *b){
  int i;
  optimizer OP;
  OP.lamb = 1.0; /* this can be updated after this function */
  OP.r = 0.5;
  OP.nx = nx;
  OP.nb = nb;
  OP.y = y;
  OP.x = x;
  OP.b = b;
  OP.f = malloc (nx*sizeof (double));
  OP.J = malloc (nx*sizeof (double*));
  for (i=0; i<nx; i++) OP.J[i] = malloc (nb*sizeof (double));
  OP.LHS1 = malloc (nb*sizeof (double*));
  OP.LHS2 = malloc (nb*sizeof (double*));
  for (i=0; i<nb; i++){
    OP.LHS1[i] = malloc (nb*sizeof (double));
    OP.LHS2[i] = malloc (nb*sizeof (double));
  }
  OP.rhs = malloc (nb*sizeof (double));
  OP.b1 = malloc (nb*sizeof (double));
  OP.b2 = malloc (nb*sizeof (double));
  return OP;
}

/* --------------------------------------------------------------------------
    Perform one step and return errors for two guesses
   -------------------------------------------------------------------------- */
void
lev_marq_step (optimizer OP, double *err1, double *err2){
  int i, j, k;
  lev_marq_function (OP.f, OP.x, OP.b, OP.nx, OP.nb);
  lev_marq_jacobian (OP.J, OP.x, OP.b, OP.nx, OP.nb);
  for (i=0; i<OP.nb; i++){
    for (j=0; j<OP.nb; j++){
      OP.LHS1[i][j] = 0.0;
      for (k=0; k<OP.nx; k++){
        OP.LHS1[i][j] += OP.J[k][i]*OP.J[k][j];
      }
      OP.LHS2[i][j] += OP.LHS1[i][j];
    }
    OP.LHS1[i][i] -= OP.lamb*OP.LHS1[i][i];
    OP.LHS2[i][i] -= OP.r*OP.lamb*OP.LHS1[i][i];
  }
  for (i=0; i<OP.nb; i++){
    OP.rhs[i] = 0.0;
    for (k=0; k<OP.nx; k++){
      OP.rhs[i] += OP.J[k][i]*(OP.y[k]-OP.f[k]);
    }
  }
  gauss_solve (OP.LHS1, OP.b1, OP.rhs, OP.nb);
  gauss_solve (OP.LHS2, OP.b2, OP.rhs, OP.nb);
  for (i=0; i<OP.nb; i++){
    OP.b1[i] += OP.b[i];
    OP.b2[i] += OP.b[i];
  }
  lev_marq_function (OP.f, OP.x, OP.b1, OP.nx, OP.nb);
  *err1 = 0.0;
  for (k=0; k<OP.nx; k++){
    *err1 += (OP.y[k]-OP.f[k])*(OP.y[k]-OP.f[k]);
  }
  lev_marq_function (OP.f, OP.x, OP.b2, OP.nx, OP.nb);
  *err2 = 0.0;
  for (k=0; k<OP.nx; k++){
    *err2 += (OP.y[k]-OP.f[k])*(OP.y[k]-OP.f[k]);
  }
}

/* --------------------------------------------------------------------------
    Optimize parameters until relative changes in the error fall below
    derr_min or the max number of steps is reached.
   -------------------------------------------------------------------------- */
void
lev_marq_optimize (optimizer OP, double *err, double derr_min,
                   int max_steps){
  int k, step;
  double err0, err1, err2, derr, *tmp;
  /* initialize error */
  lev_marq_function (OP.f, OP.x, OP.b, OP.nx, OP.nb);
  err0 = 0.0;
  for (k=0; k<OP.nx; k++){
    err0 += (OP.y[k]-OP.f[k])*(OP.y[k]-OP.f[k]);
  }
  derr = 2.0*derr_min;
  for (step = 0; step<max_steps && derr > derr_min; step++){
    lev_marq_step (OP, &err1, &err2);
    //printf ("err0 = %.5e  err1 = %.5e  err2 = %.5e\n",
    //        err0, err1, err2);
    if (err1 < err0 && err1 < err2){
      tmp = OP.b; OP.b = OP.b1; OP.b1 = tmp;
      derr = (err0-err1)/err0;
      err0 = err1;
    }
    else if (err2 < err0 && err2 < err1){
      tmp = OP.b; OP.b = OP.b2; OP.b2 = tmp;
      derr = (err0-err2)/err0;
      err0 = err2;
      OP.lamb = OP.r*OP.lamb;
    }
    else if (err2 < err1) OP.lamb = OP.r*OP.lamb;
    else OP.lamb = OP.lamb/OP.r;
  }
  *err = err0;
}
