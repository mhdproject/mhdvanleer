
#ifndef ROE_CLASS
#define ROE_CLASS
#include "global.h"
class Roe {
  ofstream outFile;
  int ii, jj;
  int d[2];

  double rr, ri, px, py, et, ke;
  double rl, ul, vl, pl, hl, al;
  double ur, vr, pr, hr, ar;
  double rho_rl, url, vrl, prl, arl, hrl;
  double kx = 0, ky = 0;
  double rho_i_2_c = 0;
  double delta_p;
  double delta_v;
  double delta_u;
  double delta_rho;

  double *lambda;
  double *cc;
  double *eigenwt;
  double *lres_state_prim;
  double *rres_state_prim;
  double *delta_w;
  double *rflux;
  double *lflux;
 public:
  int roe(const double *leftstate, const double *rightstate, double *flux, int iii, int jjj, int idir);
};
#endif
