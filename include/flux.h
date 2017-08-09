#include "global.h"

class FluxCalc {
 public:

  int flux(Array3D<zone> grid, double *fp, double *fn, int iii, int jjj, int idir, int second);
  int ctop(double *c, double *p);
  int ptoc(double *p, double *c);
};

