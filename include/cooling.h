#ifndef COOLING_FUNCTION
#define COOLING_FUNCTION
#include "global.h"
#include "tabfind.h"
#include "molcool.h"
class RadiativeCooling {
 public:
  int cooling(
      const double *zone,
      double *Lcooling,
      double dt
  );
};
#endif

