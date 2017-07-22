
#include "global.h"
#include "flux.h"
#include "cooling.h"

int update(
    Array3D<zone> newgrid,
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
    double dt,
    int second
);
