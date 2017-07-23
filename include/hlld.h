#include "global.h"
class Hlld {
 public:
  int hlld(
      const double *leftstate,
      const double *rightstate,
      double *fhlld,
      double *Res_state,
      int time_step,
      double *max_speed,
      int idir);
};

