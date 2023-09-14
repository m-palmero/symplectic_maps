# ifndef _symplectic_maps_analysis_h
# define _symplectic_maps_analysis_h

// ----------------- Dependencies ----------------------
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <time.h>
# include <omp.h>
# include "../lib/auxiliary_functions.c"
# include "../lib/gnuplot_functions.c"
# include "../lib/asprintf.c"
// --------------- Hamiltonian Maps --------------------
# include "maps/standard_map.c"
# include "maps/boozer_map.c"
# include "maps/ullmann_map.c"
# include "maps/fermi_ulam_map.c"

struct map
{
    double parameter;
    double x0_center;
    double y0_center;
    double x0_min;
    double x0_max;
    double y0_min;
    double y0_max;
    double x_space_min;
    double x_space_max; 
    double y_space_min; 
    double y_space_max;
};

typedef struct map MAP;

# endif

