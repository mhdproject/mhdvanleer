#include "global.h"
#include "vanleer.h"
#include "roe.h"
#include "riemann.h"
#include "vanleer_fvsplit.h"

class FluxCalc {
 public:

  int flux(Array3D<zone> grid, double *fp, double *fn, int iii, int jjj, int idir, int second);
};

