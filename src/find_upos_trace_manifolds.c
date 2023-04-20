// AUTOMAN DEPENDENCES
# include "automan/matrix.c"
# include "automan/map_jacobian_n.c"
# include "automan/gauss_solver.c"
# include "automan/periodic_orbit.c"

// USER CONTROLS
# include "global_parameters.c"

// AUTOMAN
# include "automan/automan.c"

// USER FUNCTIONS
# include "map_for_automan.c"

# include "dynamics.h"

int find_upos_trace_manifolds(MAP *m)
{
    char filename[T];
    manifold M;
    double x0[2];     // fixed point or member of upo
    double a = m->parameter;

    // initial guess for the periodic orbit coordinates
    x0[0] = m->x0_center;
    x0[1] = m->y0_center;

    // search orbit from the guess
	periodic_orbit(PERIOD, x0, &a);
    sprintf(filename, "upo.dat");
	save_po(filename, a, x0, PERIOD);

    // calculate manifold from x0
    //sprintf(filename, "manifold_UPO4.dat");
	//trace_manifold(x0, &a, &M);
    // save manifold modes to file
	//save_nodes(filename, M);

    return 0;
}