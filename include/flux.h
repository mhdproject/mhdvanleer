#ifndef FLUXCALC_CLASS
#define FLUXCALC_CLASS
#include "global.h"

class FluxCalc {
  int d[2];
  int kk = 0;
  int hh = 0;
  int status;

  double *leftstate;
  double *rightstate;

  double *pll;
  double *plm;
  double *plr;
  double *prr;
  double *leftprim;
  double *rightprim;
 public:

  int flux(Array3D<zone> grid, double *fp, double *fn, int iii, int jjj, int idir, int second);
  int ctop(double *c, double *p);
  int ptoc(double *p, double *c);
};
#endif

