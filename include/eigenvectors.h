#include "global.h"

class Eigen {

 public:
  int eigenvectors(const double *av_state, double **, double **, double **, double **, double, double);
  int eigenvalues(const double *sta, double *eigenval);
};
