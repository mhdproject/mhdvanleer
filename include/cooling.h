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
 public:
  int cooling(
      const double *zone,
      double *Lcooling,
      double dt
  );
};
#endif

