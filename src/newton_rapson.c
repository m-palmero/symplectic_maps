/* -----------------------------------------------------------------
    User defined functions
   ----------------------------------------------------------------- */
double newt_rap_func (double x, double *a);
double newt_rap_dfdx (double x, double *a);

/* -----------------------------------------------------------------
    Find fixed point from an initial guess to a given accuracy
   ----------------------------------------------------------------- */
void
newton_rapson (double *x, double *a, double tol, int *err){
  int step1, step2, max_steps;
  double x0, x1, f0, f1;
  x0 = *x;
  f0 = newt_rap_func (x0, a);
  max_steps = 20;
  step1 = 0; step2 = 0;
  while (fabs(f0) > tol && step1 < max_steps && step2 < max_steps){
    x1 = x0 - newt_rap_func (x0, a)/newt_rap_dfdx (x0, a);
    f1 = newt_rap_func (x1, a);
    step2 = 0;
    while (fabs(f1) > fabs(f0) && step2<max_steps){ // a safe cycle
      x1 = 0.5*(x0+x1);
      f1 = newt_rap_func (x1, a);
      step2++;
    }
    x0 = x1;
    f0 = f1;
    step1++;
    //printf ("f1 = %.10f\n",f1);
  }
  if (step1==max_steps || step2 == max_steps){
      //printf ("max. steps reached at x = %5.e\n", x0);
    *err = 1;
  }
  *x = x0;
}
