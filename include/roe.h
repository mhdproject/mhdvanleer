
#ifndef ROE_CLASS
#define ROE_CLASS
#include "global.h"
class Roe {
 public:
  int roe(const double *leftstate, const double *rightstate, double *flux, int iii, int jjj, int idir);
};
#endif
