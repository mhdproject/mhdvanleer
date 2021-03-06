#include "flux.h"
#define ROE
#undef VANLEER
#define SECOND_ORDER_SPACE
//#undef SECOND_ORDER_SPACE
int
FluxCalc::flux(Array3D<zone> oldgrid,
               double *InterfaceFlux,
               double *ResolvedState,
               int ii,
               int jj,
               int idir,
               int second) {


  pll = new double[ne];
  plm = new double[ne];
  plr = new double[ne];
  prr = new double[ne];
  leftprim = new double[ne];
  rightprim = new double[ne];

  leftstate = new double[ne];
  rightstate = new double[ne];

  if (idir == 1) {
    d[0] = 1;
    d[1] = 0;
  } else if (idir == 2) {
    d[0] = 0;
    d[1] = 1;
  }

  /* Determine the x fluxes */
  for (hh = 0; hh < ne; hh++) {
    leftstate[hh] = oldgrid[ii - d[0]][jj - d[1]][kk].array[hh];
    rightstate[hh] = oldgrid[ii][jj][kk].array[hh];
//               leftstate[hh]  = fabs(leftstate[hh]) > 1e-10?  leftstate[hh] :0; 
//               rightstate[hh]  = fabs(rightstate[hh]) > 1e-10?  rightstate[hh] :0; 
//               cout << "flux:" << leftstate[hh] << " " << rightstate[hh] << endl;
  }

  status = ctop(leftstate, leftprim);
  assert(status == 0);

  status = ctop(rightstate, rightprim);
  assert(status == 0);

//     Need a 2nd order in space correction and a non-linear flux limiter here

#ifdef SECOND_ORDER_SPACE
  if (second == 1)
//  if (1)
  {

    status = ctop(oldgrid[ii - 2 * d[0]][jj - 2 * d[1]][kk].array, pll);
    assert(status == 0);

    if (status == 1) {
      cout << "Location: " << ii << ", " << jj << endl;
      exit(0);
    }
    pll[4] = (max)(pmin, pll[4]);
    status = ctop(oldgrid[ii - d[0]][jj - d[1]][kk].array, plm);
    assert(status == 0);

    if (status == 1) {
      cout << "Location: " << ii << ", " << jj << endl;
      exit(0);
    }
    plm[4] = (max)(pmin, plm[4]);
    status = ctop(oldgrid[ii][jj][kk].array, plr);
    assert(status == 0);
    if (status == 1) {
      cout << "Location: " << ii << ", " << jj << endl;
      exit(0);
    }
    plr[4] = (max)(pmin, plr[4]);
    status = ctop(oldgrid[ii + d[0]][jj + d[1]][kk].array, prr);
    assert(status == 0);

    if (status == 1) {
      cout << "Location: " << ii << ", " << jj << endl;
      exit(0);
    }
    prr[4] = (max)(pmin, prr[4]);

    for (hh = 0; hh < ne; hh++) {
      double slope1;
      double slope2;
      double mid;
      double right;
      double left;
//               left  = oldgrid[ii - 2 * d[0]][jj - 2 * d[1]][kk].array[hh];
//               mid   = oldgrid[ii -     d[0]][jj - d[1]][kk].array[hh];
//               right = oldgrid[ii           ][jj][kk].array[hh];
      left = pll[hh];
      mid = plm[hh];
      right = plr[hh];
      slope1 = vanleer((mid - left), (right - mid));

//               left  = oldgrid[ii - d[0]][jj - d[1]][kk].array[hh];
//               mid   = oldgrid[ii       ][jj       ][kk].array[hh];
//               right = oldgrid[ii + d[0]][jj + d[1]][kk].array[hh];
      left = plm[hh];
      mid = plr[hh];
      right = prr[hh];
      slope2 = vanleer((mid - left), (right - mid));

      leftprim[hh] = plm[hh] + 0.5 * delta_x * slope1;
      rightprim[hh] = plr[hh] - 0.5 * delta_x * slope2;
    }

    // Flatten slopes here?

    status = ptoc(leftprim, leftstate);
    assert(status == 0);

    status = ptoc(rightprim, rightstate);
    assert(status == 0);
    if (leftprim[4] < 1e-10) {
      cout << "leftprim[4] < 1e-10 " << leftprim[4]
           << " " << oldgrid[ii - d[0]][jj - d[1]][kk].array[hh]
           << " " << ii << " " << jj << endl;
      if (leftprim[4] < 0) {
        cout << " Exiting on pressure, " << leftprim[4] << " < 0 " <<
             endl;
        exit(0);
      }
    }
    if (rightprim[4] < 1e-10) {
      cout << "rightprim[4] < 1e-10 " << rightprim[4]
           << " " << oldgrid[ii][jj][kk].array[hh]
           << " " << ii << " " << jj << endl;
      if (rightprim[4] < 0) {
        cout << " Exiting on pressure, " << rightprim[4] << " < 0 " <<
             endl;
        exit(0);
      }
    }

  }
#endif

#ifdef ROE
  if (leftstate[0] < 0 || rightstate[0] < 0) {
    cout << "flux:negative density in roe solver" << endl;
    cout << "Location = " << ii << "," << jj << endl;
    exit(0);
  }

  // Rotate fluxes here
  double lstate[8];
  double rstate[8];
  double iflux[8];
  double Res_state[8];

  if (idir == 1) {
    for (hh = 0; hh < ne; hh++) {
      lstate[hh] = leftstate[hh];
      rstate[hh] = rightstate[hh];
    }
  } else if (idir == 2) {
    rstate[0] = rightstate[0];
    rstate[1] = rightstate[2];
    rstate[2] = rightstate[3];
    rstate[3] = rightstate[1];
    rstate[4] = rightstate[4];
    rstate[5] = rightstate[6];
    rstate[6] = rightstate[7];
    rstate[7] = rightstate[5];

    lstate[0] = leftstate[0];
    lstate[1] = leftstate[2];
    lstate[2] = leftstate[3];
    lstate[3] = leftstate[1];
    lstate[4] = leftstate[4];
    lstate[5] = leftstate[6];
    lstate[6] = leftstate[7];
    lstate[7] = leftstate[5];
  }


//status = solver (leftstate, rightstate, InterfaceFlux, ResolvedState, timestep, &unused, idir);
  Hlld hlld;
  status = hlld.solver(leftprim, rightprim, iflux);
  assert(status == 0);

  if (status == 1) {
    cout << "location " << ii << jj << endl;
    exit(0);
  }

  if (idir == 1) {
    for (hh = 0; hh < ne; hh++) {
      InterfaceFlux[hh] = iflux[hh];
      ResolvedState[hh] = Res_state[hh];
    }
  } else if (idir == 2) {
    InterfaceFlux[0] = iflux[0];
    InterfaceFlux[2] = iflux[1];
    InterfaceFlux[3] = iflux[2];
    InterfaceFlux[1] = iflux[3];
    InterfaceFlux[4] = iflux[4];
    InterfaceFlux[6] = iflux[5];
    InterfaceFlux[7] = iflux[6];
    InterfaceFlux[5] = iflux[7];

    ResolvedState[0] = Res_state[0];
    ResolvedState[2] = Res_state[1];
    ResolvedState[3] = Res_state[2];
    ResolvedState[1] = Res_state[3];
    ResolvedState[4] = Res_state[4];
    ResolvedState[6] = Res_state[5];
    ResolvedState[7] = Res_state[6];
    ResolvedState[5] = Res_state[7];
  }

  // Unrotate fluxes

  if (status != 0) {
    cout << "Location = " << ii << "," << jj << endl;
    exit(0);
  }


  // Artificial Viscosity goes in here  ----
  //
  // LAPIDUS method divv()*q * Conserved Differences* 
  //
  // InterfaceFlux[0]=InterfaceFlux[0]+max(0,)


#ifdef LAPIDUS_ARTIFICIAL_VISCOSITY
  double    u1 = oldgrid[ii + 1][jj    ][0] _MOMX/oldgrid[ii + 1][jj    ][0] _MASS;
  double    u2 = oldgrid[ii - 1][jj    ][0] _MOMX/oldgrid[ii - 1][jj    ][0] _MASS;
  double    v1 = oldgrid[ii    ][jj + 1][0] _MOMY/oldgrid[ii    ][jj + 1][0] _MASS;
  double    v2 = oldgrid[ii    ][jj - 1][0] _MOMY/oldgrid[ii    ][jj - 1][0] _MASS;
  double    divv = 0.5 * (u1 - u2 +v1 -v2);


   for (hh=0 ; hh<ne ; hh++)
   {
  double    delu = 0;
   delu=(oldgrid[ii][jj][0].array[hh]-oldgrid[ii -d[0]][jj -d[1] ][0].array[hh] );
  InterfaceFlux[hh]=InterfaceFlux[hh]+ dtodx*0.1*i(min)(0.0,divv)*delu;
   }
#endif  /* LAPIDUS_ARTIFICIAL_VISCOSITY */

#endif /* ROE */

  delete[] leftstate;
  delete[] rightstate;
  delete[] pll;
  delete[] plm;
  delete[] plr;
  delete[] prr;
  delete[] leftprim;
  delete[] rightprim;

  return 0;

}

int
FluxCalc::ptoc(double *p, double *c) {
  double rho;
  double velx, vely, velz, velx2, vely2, velz2, ke, pressure;
  double bx, by, bz, bsquared;
  double gammam1 = gammag - 1;
  double gammam1i = 1 / gammam1;
  rho = p[0];
  velx = p[1];
  vely = p[2];
  velz = p[3];
  pressure = p[4];
  bx = p[5];
  by = p[6];
  bz = p[7];

  velx2 = velx * velx;
  vely2 = vely * vely;
  velz2 = velz * velz;
  ke = 0.5 * rho * (velx2 + vely2 + velz2);
  bsquared = bx * bx + by * by + bz * bz;
  c[0] = rho;
  c[1] = rho * velx;
  c[2] = rho * vely;
  c[3] = rho * velz;
  c[4] = ke + pressure * gammam1i + 0.5 * bsquared;
  c[5] = bx;
  c[6] = by;
  c[7] = bz;
  /*
     cout << "p " << p[0] 
     << " " << p[1] 
     << " " << p[2] 
     << " " << p[3] 
     << endl;
     cout << "c " << c[0] 
     << " " << c[1] 
     << " " << c[2] 
     << " " << c[3] 
     << endl;
   */

  if (c[4] <= 1e-20) {
    cout << "p " << p[0]
         << " " << p[1] << " " << p[2] << " " << p[3] << endl;
    cout << "c " << c[0]
         << " " << c[1] << " " << c[2] << " " << c[3] << endl;
    cout << "ptoc: Mein Gott!" << endl;
//              exit (0);
  }

  if (vely < -400) {
    cout << setiosflags(ios::fixed);
    cout << "ptoc: p " << p[0] << " " << p[1] << " " << p[2] << " " << p[3]
         << endl;
    cout << "ptoc: c " << c[0] << " " << c[1] << " " << c[2] << " " << c[3]
         << endl;
  }

  return 0;
}

int
FluxCalc::ctop(double *c, double *p) {
  double rho, rhoi;
  double px;
  double py, pz, et, velx, vely, velz, velx2, vely2, velz2, ke, pressure;
  double bx, by, bz, bsquared;
  double gammam1 = gammag - 1;
  rho = c[0];
  px = c[1];
  py = c[2];
  pz = c[3];
  et = c[4];
  bx = c[5];
  by = c[6];
  bz = c[7];
  rhoi = 1.0 / rho;
  velx = px * rhoi;
  vely = py * rhoi;
  velz = pz * rhoi;
  velx2 = velx * velx;
  vely2 = vely * vely;
  velz2 = velz * velz;
  ke = 0.5 * rho * (velx2 + vely2 + velz2);
  bsquared = bx * bx + by * by + bz * bz;
  pressure = et - ke - 0.5 * bsquared;
  pressure = pressure * (gammam1);
  p[0] = rho;
  p[1] = velx;
  p[2] = vely;
  p[3] = velz;
  p[4] = pressure;
  p[5] = bx;
  p[6] = by;
  p[7] = bz;

  if (pressure < 1e-20) {
    cout << "p " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] <<
         endl;
    cout << "c " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] <<
         endl;
    cout << "ctop: pressure low" << endl;
    if (pressure < 0) {
      cout << "pressure <0 " << pressure << endl;
      return 1;
    }
  }
  if (vely < -400) {
    cout << setiosflags(ios::fixed);
    cout << "Vel < 400" << endl;
    cout << "p " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] <<
         endl;
    cout << "c " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] <<
         endl;
  }
  return 0;
}
