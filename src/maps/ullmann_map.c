/* ----------------------------------------------------------
    PARAMETERS
    a[0] = C // perturbation strength
    a[1] = m // perturbation mode
    a[2] = y_a = 1-a/b // dimensionless plasma radius
    a[3] = q_a // plasma edge safety factor
    a[4] = a1 // toroidicity
    a[5] = gamma // iota parameter
    a[6] = y_star // used to carry auxiliary values
    a[7] = x_star
   ---------------------------------------------------------- */

double heaviside (double x)
{
  if (x>=0.) return 1.;
  else return 0.;
}

double iota (double y, double *a)
{
  double u;
  u = (1.0 - y)/(1.0 - a[2]); u *= u;
  return (1.-pow(1.-u, a[5]+1.)*heaviside(y-a[2]))/(a[3]*u);
}

double newt_rap_func (double y, double *a)
{
  return y-a[6]-(a[1]*a[0]/(a[1]-1.))*pow(1.-y, a[1]-1.)*sin (a[1]*a[7]);
}

double newt_rap_dfdx (double y, double *a)
{
  return 1.+a[1]*a[0]*pow(1.-y,a[1]-2.)*sin(a[1]*a[7]);
}

void newton_rapson (double *x, double *a, double tol, int *err)
{
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

/* ----------------------------------------------------------
    The Ulmann map
   ---------------------------------------------------------- */
void ullmann_map_eqs(double *x0, double *x1, double *parameter)
{
    int mode, err;
    double r_a, r_b, q_a, dB_B;
    double *a;
    double *xs, ****x;

    double *eps = (double*) parameter;
    
    /* Physical parameters */
    mode = 7;     
    q_a = 5.0;    
    r_a = 0.18;   
    r_b = 0.22;  
    //eps = 0.1/100;
    dB_B = *parameter; 


    /* Ullmann map parameters */
    a = malloc (8*sizeof(double));
    a[0] = 2.*M_PI*pow(r_b/r_a,mode-2)*dB_B/q_a; // pert. strenght
    a[1] = 1.0*mode;      // perturbation mode
    a[2] = 1.0 - r_a/r_b; // plasma edge loc. relative to limiter
    a[3] = q_a;
    a[4] =-0.04; // toroidicity
    a[5] = 3.0;   // gamma, pol. field parameter
    a[6] = 0.0;   // y_star // used to carry auxiliary values
    a[7] = 0.0;   // x_star
    
    err = 0;
    a[6] = 1.-(1.-x0[1])/(1.-a[4]*sin(x0[0]));
    a[7] = x0[0] + 2.*M_PI*iota(a[6], a) + a[4]*cos(x0[0]);
    // initialize y_n+1 = y_star
    x1[1] = a[6] + a[0]*pow(1.-a[6], a[1]-1)*sin(a[1]*a[7]); 
    newton_rapson (x1+1, a, 1e-12, &err); // find y_n+1 = x1[1] = x1+1
    if (err == 0)
    {
	    x1[0] = a[7] - a[0]*pow(1.-x1[1],a[1]-2.)*cos(a[1]*a[7]);
	    x1[0] = fmod(x1[0],2.*M_PI);
    }
    else 
    { // return nans
	    //printf ("orbit exit domain\n");
	    x1[0] = NAN;
	    x1[1] = NAN;
    }
    // x0[0] = x1[0];
    // x0[1] = x1[1];
    free(a);
}

