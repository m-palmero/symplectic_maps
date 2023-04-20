/* extermal dependences */
// gauss_solver.c

/* --------------------------------------------------------------------------
    A structure that encapsulates all the required information
   -------------------------------------------------------------------------- */
typedef struct {
  double lamb; /* damping parameter */
  double r;   /* damping adjustment parameter */
  double **J; /* the jacobian of the levenberg-marquardt function */
  double *f;  /* the vector of modeled values */
  double *y;  /* the set of known dependent values */
  double *x;  /* the set of known independent variables */
  double *b;  /* the set of optimization parameters */
  int nx;     /* the number of model points */
  int nb;     /* the number of optimization parameters */
  /* auxiliary arrays */
  double **LHS1; /* the lhs matrix */
  double **LHS2; /* the lhs matrix */
  double *rhs; /* the rhs vector */
  double *b1; /* param variation for lamb */
  double *b2; /* param variation for r*lamb */
} optimizer;
