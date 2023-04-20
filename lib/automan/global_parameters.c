/* ----------------------------------------------------
   Global variables
   ---------------------------------------------------- */
/* UPO period */ int
PERIOD = 31;

/* Manifold type: U, S */ char
MAN_TYPE = 'U';

/* manifold branch: L, R */ char
BRANCH = 'L';

/* maximum chord_length */ double
MAX_CHD = 1e-2;

/* minimum chord_length */ double
MIN_CHD = 1e-3;

/* maximum chord angle in deg. */ double
MAX_ANG = 5.;

/* number of segments  */ int
N_SEG = 10;

/* dist from first segment to saddle */ double
DS = 1e-7;

/* space added during memory reallocation */ int
MEM_BLOCK = 100;

/* a global map call counter */ int
MAP_CALL_COUNT = 0;

/* computational box */
double X_MIN = -10.0;
double X_MAX = 10.0;
double Y_MIN = -10.0;
double Y_MAX = 10.0;

