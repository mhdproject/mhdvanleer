
#ifndef RIEMANN_CLASS
#define RIEMANN_CLASS
#include "global.h"
class Riemann {
 public:
  int solver(double *leftstate, double *rightstate, double *roeflux, double *res_state, int idir);
};
#endif
