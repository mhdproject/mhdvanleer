#include "global.h"
#include "vanleer.h"
#include "roe.h"
#include "riemann.h"
#include "vanleer_fvsplit.h"

int flux(Array3D<zone> grid,
         double *fp,
         double *fn,
         double dtodx,
         int iii,
         int jjj,
         int timestep,
         int idir,
         int second);

