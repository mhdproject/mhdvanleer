#include "global.h"
int initialise(char *filename, Array3D<zone> grid, int *, double *);

class InitialJet {

 public:
  int setup(char *filename, Array3D<zone> grid, int *, double *);
};

class InitialBlast {

 public:
  int setup(char *filename, Array3D<zone> grid, int *, double *);
};

class InitFactory {
 public:
  int make_init(string probtype, int argc, Array3D<zone> grid, int maxstep, double cfl, char **argv);

};
