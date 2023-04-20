/* --------------------------------------------------------------------------
    Solve Ax = b for x using Gaussian elimination
   -------------------------------------------------------------------------- */
void
gauss_solve (double** A, double* x, double* b, int dim){
  int i, j, k;
  double temp;

  bool back_sub = true;

  double** Ab; /* augmented matrix */
  Ab = malloc (dim*sizeof(double*));

  for (i=0; i<dim; i++){
    Ab[i] = malloc ((dim+1)*sizeof(double));
    for (j=0; j<dim; j++) Ab[i][j] = A[i][j];
    Ab[i][dim] = b[i];
  }
  for (k=0; k<dim; k++){ /* move in the pivots */
    if (Ab[k][k] == 0.0){ /* if zero swaps with a nonzero */
      i = k+1;
      while (Ab[i][k] == 0 && i<dim) i++;
      if (i == dim){
	       back_sub = false;
	        break;
      }
      for (j=0; j<dim+1; j++){
	       temp = Ab[k][j];
	        Ab[k][j] = Ab[i][j];
	        Ab[i][j] = temp;
      }
    }
    for (j=dim; j>=k; j--) Ab[k][j] /= Ab[k][k];
    for (i=k+1; i<dim; i++){
      if (Ab[i][k] != 0.0){
	       for (j=dim; j>=k; j--){
	          Ab[i][j] = Ab[i][j] - Ab[i][k]*Ab[k][j];
	       }
      }
    }
  }
  /* performs the back_substitution */
  if (back_sub){
    for (i=dim-1; i>=0; i--){
      temp = Ab[i][dim];
      for (j = i+1; j<dim; j++)	temp -= Ab[i][j]*x[j];
      x[i] = temp;
    }
  }
} /* end of gauss_solver */
