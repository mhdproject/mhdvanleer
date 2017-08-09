#ifndef COOLING_CLASS
#define COOLING_CLASS
#include "global.h"
#include "tabfind.h"
#include "molcool.h"
class RadiativeCooling {
  double gammam1 = gammag - 1;
  double rhoi, px, py, et, ke;
  double rho, velx, vely, pressure;
  double velx2;
  double vely2;
  double temperature;

  double ki = 24296.3696;
  double mp = 1.67262158;
  double mpi = 1.0 / mp;
  double nt = 0;

  double rate = 0, nt2, de, y, subdt;
  double eloss;
  double molrate = 0;
  double lowest_temperature, subttot;
  double e_init;
  int firststep;
  int counter;
  double nH2;
  double nh;
  double chi;
  double real_temp;
  int rc = 0;

  double atomic_temp_tab[50] = {
      4.04, 4.08, 4.12, 4.17, 4.21, 4.25, 4.3, 4.34, 4.38, 4.44, 4.48, 4.53,
      4.58, 4.62, 4.67, 4.71, 4.77, 4.82, 4.88, 4.93, 4.98, 5.04, 5.09, 5.14,
      5.2, 5.25, 5.3, 5.36, 5.41, 5.46, 5.51, 5.57, 5.62, 5.67, 5.73, 5.78,
      5.83, 5.88, 5.94, 5.99, 6.04, 6.09, 6.14, 6.2, 6.25, 6.3, 6.35, 6.4,
      6.45, 6.5
  };
 public:
  int cooling(
      const double *zone,
      double *Lcooling,
      double dt
  );
};
#endif

