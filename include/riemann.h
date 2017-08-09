
#ifndef RIEMANN_CLASS
#define RIEMANN_CLASS
#include "global.h"
class Riemann {
  ofstream outFile;
  char outfilename[50] = "output/solver.log";

  int ii = 0;
  int jj = 0;
  int kk = 0;

  double av_state[8];
  double cslow = 0;
  double cslow2 = 0;
  double calfven = 0;
  double calfven2 = 0;
  double cfast = 0;
  double cfast2 = 0;
  double bsquared = 0;

  double rho_rl = 0;
  double u_rl = 0;
  double v_rl = 0;
  double w_rl = 0;
  double p_rl = 0;
  double bu_rl = 0;
  double bv_rl = 0;
  double bw_rl = 0;

  double lambda[7];
  double lstate[7];
  double rstate[7];
  double lrsp[7];
  double rrsp[7];
  double eigenwt[7];

  double rho = 0;
  double mass = 0;
  double rhoi = 0;
  double u = 0;
  double v = 0;
  double w = 0;
  double bu = 0;
  double bv = 0;
  double bw = 0;
  double p = 0;
  double csound2 = 0;
  double term = 0;
  double a_star2 = 0;
  double kinetic = 0;
  double p_magnetic = 0;
  double internal_energy = 0;
  double energy = 0;
  double v2 = 0;
  double vdotb = 0;

  double cc[7];
  double dvy = 0, dvz = 0;
 public:
  int solver(double *leftstate, double *rightstate, double *roeflux, double *res_state, int idir);
};
#endif
