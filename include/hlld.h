#ifndef  HLLD_CLASS
#define  HLLD_CLASS
#include "global.h"
class Hlld {
 public:
  int solver(const double *leftstate, const double *rightstate, double *fhlld);
};
#endif

