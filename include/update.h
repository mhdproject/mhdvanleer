
#include "global.h"
#include "flux.h"
#include "cooling.h"

class Updater {

 public:
  int update(Array3D<zone> newgrid,
             Array3D<zone> oldgrid,
             Array3D<zone> fx,
             Array3D<zone> fy,
             Array3D<zone> xResState,
             Array3D<zone> yResState,
             double del,
             int ii,
             int jj,
             int timestep,
             Array3D<zone> fluxgrid,
             double dt);

  Array3D<zone> &check_negative_pressure(Array3D<zone> &NewGrid,
                                         const Array3D<zone> &oldg,
                                         int ii,
                                         int jj,
                                         int hh,
                                         int kk,
                                         double ke,
                                         const double *fx1,
                                         const double *fx2,
                                         const double *fy1,
                                         const double *fy2,
                                         double b2,
                                         double pressure) const;
};


