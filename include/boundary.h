#ifndef BOUNDARY_CLASS
#define BOUNDARY_CLASS
#include "global.h"
class BoundaryCondition {
 public:
  int boundary(Array3D<zone> grid, int inject_jet);
};
#endif

