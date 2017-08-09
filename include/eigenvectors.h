#ifndef  EIGEN_CLASS
#define  EIGEN_CLASS
#include "global.h"
#include <boost/math/constants/constants.hpp>

class Eigen {
  int ii = 0;
  int jj = 0;
  double rho = 0;
  double rhoi = 0;
  double u = 0;
  double v = 0;
  double w = 0;
  double BBv = 0;
  double BBw = 0;
  double bu = 0;
  double bv = 0;
  double bw = 0;
  double p = 0;
  double cslow = 0;
  double calfven = 0;
  double calfven2 = 0;
  double cfast = 0;
  double cslow2 = 0;
  double cfast2 = 0;
  double csound2 = 0;
  double csound = 0;
  double term = 0;
  double a_star2 = 0;
  double bsquared = 0;
  double vv2 = 0;

  double betay = 0;
  double betaz = 0;
  double alphaf = 0;
  double alphas = 0;
  double bperp = 0;
  double icsound2 = 0;
  double icsound22 = 0;
  double phi = 0;
  double deltas = 0;
  double deltaf = 0;
  double sqrt_rho = 0;
  double sqrt_rhoi = 0;
  double pie;

 public:
  int eigenvectors(const double *av_state, double **, double **, double **, double **, double, double);
  int eigenvalues(const double *sta, double *eigenval);
};
#endif
