/* ============== LIBRARY REQUIREMENTS ===============
   # include "automan_tools/matrix.c"
   # include "automan_tools/map_jacobian_n.c"
   # include "automan_tools/params_sample.c"
   ---------------------------------------------------- */

/* ============== USER REQUIREMENTS ==================
   User must provide a map, its inverse, jacobian and
   a method to calculate the vector from x0 to x1
   ---------------------------------------------------- */
int map_fw (double *x0, double *x1, void *pars);    // forward map
int map_bw (double *x0, double *x1, void *pars);    // inverse map
void jacobian (double **J, double *x, void *pars);  // J_mapfw
void vec_dif (double *x01, double *x1, double *x0); // x1-x0

/* ===================== USAGE ========================
   Set manifold parameters in manifold_params.c and
   include in your program. Define the forward and back
   wards mappings. Declare a manifold and call:

   trace_manifold (double *xf, void *pars, manifold *M)

   providing one member of the periodic orbit and the
   relevant mapping parameters. Save manifold to file
   using the function:

   save_nodes (char *filename, manifold M)
   ---------------------------------------------------- */


/* =========== FUNCTIONS AND STRUCTURES =============== */

/* ----------------------------------------------------
   A segment contains nodes, navigators, lengths,
   inter-secant angles, etc.
   ---------------------------------------------------- */
typedef struct {
    int space;   /* available space in arrays */
    int np;      /* number of nodes */
    int *nx;     /* the forward navigator */
    int *pv;     /* the reverse navigator */
    double *l;   /* chord lengths */
    double *cs;  /* cosine of intersecant angles */
    double **x;  /* the np nodes of the curve */
} segment;

/* ----------------------------------------------------
   The manifold is just a collection of segments
   ---------------------------------------------------- */
typedef struct {
    int ns;      /* the number of segments in the manifold */
    void *pars;  /* the system parameters */
    segment *P;  /* An array of segments */
} manifold;

/* ----------------------------------------------------
   Segment allocators
   ---------------------------------------------------- */
void
allocate_segment (int block, segment *P){
  int i;
  P->space = block;
  P->nx = malloc ((P->space)*sizeof(int));
  P->pv = malloc ((P->space)*sizeof(int));
  P->l = malloc ((P->space)*sizeof(double));
  P->cs = malloc ((P->space)*sizeof(double));
  P->x = malloc ((P->space)*sizeof(double*));
  for (i=0; i<P->space; i++){
    P->x[i] = malloc (2*sizeof(double));
  }
}

void
reallocate_segment(int block, segment *P){
  int i, space0;
  space0 = P->space;
  P->space += block;
  P->nx = realloc (P->nx, (P->space)*sizeof(int));
  P->pv = realloc (P->pv, (P->space)*sizeof(int));
  P->l = realloc (P->l, (P->space)*sizeof(double));
  P->cs = realloc (P->cs, (P->space)*sizeof(double));
  P->x = realloc (P->x, (P->space)*sizeof(double*));
  for (i=space0; i<P->space; i++){
    P->x[i] = malloc (2*sizeof(double));
  }
}

/* ----------------------------------------------------
   allocate 'block' segments in manifold
   ---------------------------------------------------- */
void
allocate_manifold (int block, manifold *M){
  M->ns = block;
  M->P = malloc (M->ns*sizeof (segment));
}

/* ----------------------------------------------------
   Main function
   ---------------------------------------------------- */
void trace_manifold (double *xf, void *pars, manifold *M);

/* ----------------------------------------------------
   Auxiliary functions
   ---------------------------------------------------- */
void first_segment (double *xf, segment *P0, void *pars);
void map_segment (int i, manifold *M);
void refine_segment (int n, manifold *M);
void new_node_after (int i, int n, manifold *M);
void map_node (int i, int n1, int n2, manifold *M);
void map (double *x0, double *x1, void *pars);
double vec_norm (double *x);
bool in_box (double *x);

/* ----------------------------------------------------
   Calculate the unstable manifold with N_SEG segments
   along the branch specified by the user and based in
   the fixed point xf.
   ---------------------------------------------------- */
void
trace_manifold (double *xf, void *pars, manifold *M){
  int i;
  
  allocate_manifold (N_SEG*PERIOD, M);
  M->pars = pars;
  /* calculate first segment */
  first_segment (xf, &(M->P[0]), M->pars);
  /* advect this segment to all branches */
  for (i=0; i<PERIOD-1; i++){
      printf ("seeding %d branch\n", i+1);
      map_segment (i, M);
  }
  /* advect all branches */
  for (i=PERIOD; i<N_SEG*PERIOD; i++){
      printf ("segment %d of %d\n", i, N_SEG*PERIOD-1);
      map_segment (i-1, M);
      // printf ("refining %d\n", i);
      refine_segment (i, M);
  }
}

/* ----------------------------------------------------
   Calculate the first segment from the fixed point a-
   long the unstable direction in the branch specified
   in the .conf file.
   ---------------------------------------------------- */
void
first_segment (double *xf, segment *P0, void *pars){
    double dx, *x0, *x1, *x, *v, **V, **J;
    
    x0 = malloc (2*sizeof(double));
    x1 = malloc (2*sizeof(double));
    
    dx = DS;
    /* the n period jacobian */
    J = matrix_alloc (2, 2);
    V = matrix_alloc (2, 2);
    v = vector_alloc (2);
    /* calculate the jacobian of T^n */
    jacobian_n (J, PERIOD, xf, pars);
    /* solve the eigenvalue problem */
    eigensystem (J, V, v);
    /* create initial condition along largest 
       eigendirection */
    if (BRANCH == 'R' && MAN_TYPE == 'U'){
	x0[0] = xf[0] + dx*V[1][0];
	x0[1] = xf[1] + dx*V[1][1];
    }
    else if (BRANCH == 'L' && MAN_TYPE == 'U'){
	x0[0] = xf[0] - dx*V[1][0];
	x0[1] = xf[1] - dx*V[1][1];
    }
    else if (BRANCH == 'R' && MAN_TYPE == 'S'){
	x0[0] = xf[0] + dx*V[0][0];
	x0[1] = xf[1] + dx*V[0][1];
    }
    else if (BRANCH == 'L' && MAN_TYPE == 'S'){
	x0[0] = xf[0] - dx*V[0][0];
	x0[1] = xf[1] - dx*V[0][1];
    }
    else printf ("Unknown branch or stability\n");
    /* find end-point of first segment */
    map_n (PERIOD, x0, x1, pars);
    /* initialize first segment as a straight segment */
    P0->np = 3;
    allocate_segment (P0->np, P0);
    
    P0->x[0][0] = x0[0]; P0->x[0][1] = x0[1];
    P0->x[2][0] = 0.5*(x0[0]+x1[0]); P0->x[2][1] = 0.5*(x0[1]+x1[1]);
    P0->x[1][0] = x1[0]; P0->x[1][1] = x1[1];
    
    P0->nx[0] = 2; P0->pv[0] = -1;
    P0->nx[2] = 1; P0->pv[2] = 0;
    P0->nx[1] = -1; P0->pv[1] = 2;
    
    P0->cs[0] = 1.0;
    P0->cs[2] = 1.0;
    P0->cs[1] = 0.0;
    
    P0->l[0] = 0.5*sqrt((x1[0]-x0[0])*(x1[0]-x0[0])+
			(x1[1]-x0[1])*(x1[1]-x0[1]));
    P0->l[2] = P0->l[0];
    P0->l[1] = 0.0;
    
    free (x0); free (x1);
    free (v); free (V); free (J);
}


/* ----------------------------------------------------
   Map segment P0 to P1 once
   ---------------------------------------------------- */
void
map_segment (int n, manifold *M){
    int i;
    double ll, lr, *dxl, *dxr, *tmp;
    segment *P0, *P1, *P;
    P0 = &(M->P[n]);
    P1 = &(M->P[n+1]);
    
    dxl = malloc (2*sizeof(double));
    dxr = malloc (2*sizeof(double));
    allocate_segment (P0->space, P1);
    P1->np = P0->np;
    for (i=0; i<P1->np; i++){
	/* copy navigators */
	P1->pv[i] = P0->pv[i];
	P1->nx[i] = P0->nx[i];
	/* set mapped point to nan if is out of box */
	/* if (in_box (P0->x[i])) map (P0->x[i], P1->x[i], M->pars); */
	/* else {P1->x[i][0] = NAN; P1->x[i][1] = NAN;} */
	map (P0->x[i], P1->x[i], M->pars);
	if (!(in_box (P1->x[i]))){
	    P1->x[i][0] = NAN;
	    P1->x[i][1] = NAN;
	}
    }
    if (n < PERIOD-1){/* use free boundaries in the first segment */
	vec_dif (dxl, P1->x[P1->nx[0]], P1->x[0]);
    }
    else {
	P = &(M->P[n-PERIOD+1]);
	vec_dif (dxl, P1->x[0], P->x[P->pv[1]]);
    }
    ll = vec_norm (dxl);
    for (i=0; i!=1; i=P1->nx[i]){
	vec_dif (dxr, P1->x[P1->nx[i]], P1->x[i]);
	lr = vec_norm (dxr);
	P1->l[i] = lr;
	P1->cs[i] = (dxr[0]*dxl[0]+dxr[1]*dxl[1])/(lr*ll);
	/* recycle values */
	tmp = dxl; dxl = dxr; dxr = tmp;
	ll = lr;
    }
    free (dxl);
    free (dxr);
}

/* ----------------------------------------------------
   Refine segment n to match the ressolution criteria
   ---------------------------------------------------- */
void
refine_segment (int n, manifold *M){
    bool new_left, new_right;
    int i;
    double ll, cos_thc;
    segment *P0, *P1, *P2; /* local references */
    P0 = &(M->P[0]);
    P1 = &(M->P[n-PERIOD]);
    P2 = &(M->P[n]);
    
    cos_thc = cos (MAX_ANG*M_PI/180.);
    
    i = 0; /* node index */
    while (i != 1){
	/* check if chords are large */
	if (i==0) ll = P1->l[P1->pv[1]];
	else ll = P2->l[P2->pv[i]];
	new_left = (ll > MAX_CHD);
	new_right = (P2->l[i] > MAX_CHD);
	/* check if angle is large and segs are above MIN_CHD */
	if (P2->cs[i] < cos_thc && P2->l[i] > MIN_CHD
	    && ll > MIN_CHD){
	    // printf ("cs[%d] = %.5e\n", i, P2->cs[i]);
	    if (P2->l[i] > ll) new_right = true;	    
	    else  new_left = true;
	}
	/* refine if required */
	if (new_right) {
	    // printf ("lr = %.5e\n", P2->l[i]);
	    // printf ("new right at i = %d\n", i);
	    new_node_after (i, n, M);
	}
	if (new_left){
	    // printf ("ll = %.5e\n", ll);
	    // printf ("new left at i = %d\n", i);
	    if (i==0){
		new_node_after (P1->pv[1], n-PERIOD, M);
		map_node (P1->pv[1], n-PERIOD, n, M);
	    }
	    else new_node_after (P2->pv[i], n, M);
	}
	/* go to new location */
	if (new_left && i!=0) i = P2->pv[P2->pv[i]];
	else if (!new_right && i!=1) i = P2->nx[i];
    }
}

/* ----------------------------------------------------
   create new node after node i in segment P[n]
   ---------------------------------------------------- */
void
new_node_after (int i, int n, manifold *M){
    segment *P0, *P1;
    P0 = &(M->P[0]);
    P1 = &(M->P[n]);
    /* add space when required */
    if (P0->np == P0->space) reallocate_segment (MEM_BLOCK, P0);
    if (P1->np == P1->space) reallocate_segment (MEM_BLOCK, P1);
    /* insert new point in first segment */
    P0->x[P0->np][0] = 0.5*(P0->x[i][0] + P0->x[P0->nx[i]][0]);
    P0->x[P0->np][1] = 0.5*(P0->x[i][1] + P0->x[P0->nx[i]][1]);
    P0->l[P0->np] = 0.5*P0->l[i]; P0->l[i] = 0.5*P0->l[i];
    P0->cs[P0->np] = 1.0;
    
    P0->nx[P0->np] = P0->nx[i];  P0->pv[P0->np] = i;
    P0->pv[P0->nx[i]] = P0->np;  P0->nx[i] = P0->np;
    P0->np++;
    /* map to new segment */
    map_node (P0->nx[i], 0, n, M);
}

/* ----------------------------------------------------
   Map node from one arbitrary position in segment n1
   to the corresponding location in n2 and calculate 
   curve parameters.
   ---------------------------------------------------- */
void
map_node (int i, int n1, int n2, manifold *M){
    int j, k;
    double *dxl, *dxr, *tmp, ll, lr;
    segment *P0, *P1, *P2;
    
    dxl = malloc (2*sizeof(double));
    dxr = malloc (2*sizeof(double));
    
    P0 = &(M->P[n1]);
    P1 = &(M->P[n2-PERIOD]);
    P2 = &(M->P[n2]);
    /* add space when required */
    if (P2->np == P2->space) reallocate_segment (MEM_BLOCK, P2);
    map_n (n2-n1, P0->x[i], P2->x[P2->np], M->pars);
    /* adjust navigators */
    P2->pv[P2->np] = P0->pv[i]; P2->nx[P2->np] = P0->nx[i];
    P2->pv[P0->nx[i]] = P2->np; P2->nx[P0->pv[i]] = P2->np;
    /* calculate parameters */
    j = P2->pv[P2->np];
    if (j == 0) vec_dif (dxl, P2->x[j], P1->x[P1->pv[1]]);
    else vec_dif (dxl, P2->x[j], P2->x[P2->pv[j]]);
    ll = vec_norm (dxl);
    for (k=0; k<3; k++){
	vec_dif (dxr, P2->x[P2->nx[j]], P2->x[j]);
	lr = vec_norm (dxr);
	P2->l[j] = lr;
	P2->cs[j] = (dxl[0]*dxr[0] + dxl[1]*dxr[1])/(ll*lr);
	/* recycle values */
	tmp = dxl; dxl = dxr; dxr = tmp;
	ll = lr;
	j = P2->nx[j];
	if (j==1) break;
    }
    P2->np++;
    free (dxl);
    free (dxr);
}

/* ----------------------------------------------------
   The map function based in the user definitions
   ---------------------------------------------------- */
void
map (double *x0, double *x1, void *pars){
    if (MAN_TYPE == 'U') map_fw (x0, x1, pars);
    else if (MAN_TYPE == 'S') map_bw (x0, x1, pars);
    else printf ("unknown manifold type\n");
}

/* ----------------------------------------------------
   Euclidian norm of vector
   ---------------------------------------------------- */
double
vec_norm (double *x){
    return sqrt (x[0]*x[0]+x[1]*x[1]);
}

/* ----------------------------------------------------
   Export nodes to an external file
   ---------------------------------------------------- */
void
save_nodes (char *filename, manifold M){
    int i, j, k, np;
    segment P;
    FILE *f0 = fopen (filename, "w");
    for (i=0; i<PERIOD; i++){
			for (j=0; j<N_SEG; j++){
					P = M.P[i+j*PERIOD];
					np = P.np;
					for (k=0; k!=1; k=P.nx[k]){
				fprintf (f0, "%.16e  %.16e\n", P.x[k][0], P.x[k][1]);
				// for connected spaces
				if (fabs(P.x[k][0]-P.x[P.nx[k]][0])>M_PI){
						fprintf (f0, "\n");
				}
					}
			}
			fprintf (f0, "\n\n");
    }
}

/* ----------------------------------------------------
   Determines if point is inside computational box
   ---------------------------------------------------- */
bool
in_box (double *x){
    if (x[0] < X_MAX && x[0] > X_MIN &&
	x[1] < Y_MAX && x[1] > Y_MIN) return true;
    else return false;
}

/* ----------------------------------------------------
   Export periodic orbit or fixed point to an external
   file
   ---------------------------------------------------- */

void save_po(char *filename, double a, double *x0, int period)
{
	FILE *out = fopen(filename, "w");
	double x1[2];
	for (int i = 0; i < period; i++)
	{
		map_fw(x0, x1, &a);
		fprintf(out, "%1.15e %1.15e\n", x1[0], x1[1]);
		copy(x0, x1, 2);
	}
	fclose(out);
}

void save_fp(char *filename, double *fp)
{
	FILE *out = fopen(filename, "w");
	fprintf(out, "%1.15e %1.15e\n", fp[0], fp[1]);
	fclose(out);
}