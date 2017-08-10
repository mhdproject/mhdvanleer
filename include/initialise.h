#ifndef INIT_CLASS
#define INIT_CLASS
#include "global.h"
int initialise(const char *filename, Array3D<zone> grid);

class InitialJet {

 public:
  int setup(const char *filename, Array3D<zone> grid);
};

class InitialBlast {

 public:
  int setup(const char *filename, Array3D<zone> grid);
};

class InitFactory {
 public:
  int make_init(string probtype, int argc, Array3D<zone> grid, int maxstep, char **argv);

};
#endif
