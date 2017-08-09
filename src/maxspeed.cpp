#include "maxspeed.h"
int
MaxSpeed::maxspeed(Array3D<zone> grid, double *maxspeed) {
  int ii = 0;
  int jj = 0;
  int kk = 0;

  double gammam1 = gammag - 1;



  /* Going through every cell in the grid, work out the sound speed and take the
   * max of the absolute value of the vx and vy. The sum constitutes the fastest wave speed
   * in the x or y directions and is passed back to maxspeed in the main program
   * and used to work out the time step */

  for (ii = 0; ii < nx - 1; ii++) {
#ifdef TWODIM
    for (jj = 0; jj < ny - 1; jj++)
#endif /* TWODIM */
    {
      double px = 0, py = 0, pz = 0, et = 0, ke = 0;
      double rl = 0, ul = 0, vl = 0, wl, pl = 0;
      double bv, bw, bsquared;
      double calfven2 = 0;
      double cfast = 0;
      double cfast2 = 0;
      double csound2 = 0;
      double term = 0;
      double a_star2 = 0;
      double vv2 = 0;
      double speed1 = 0;
      double speed2 = 0;
      double speed3 = 0;
      rl = grid[ii][jj][kk]_MASS;
      px = grid[ii][jj][kk]_MOMX;
      py = grid[ii][jj][kk]_MOMY;
      pz = grid[ii][jj][kk]_MOMZ;
      et = grid[ii][jj][kk]_ENER;

      double bu;
      bu = grid[ii][jj][kk]_B_X;
      bv = grid[ii][jj][kk]_B_Y;
      bw = grid[ii][jj][kk]_B_Z;

      double ri = 0;
      ri = 1.0 / rl;
      ul = px * ri;
      vl = py * ri;
      wl = pz * ri;
      vv2 = (ul * ul + vl * vl + wl * wl);
      ke = 0.5 * rl * vv2;
      bsquared = (bu * bu + bv * bv + bw * bw);
      pl = et - ke - 0.5 * bsquared;
      pl = pl * (gammam1);
      bsquared = bsquared * ri;

      calfven2 = bu * bu * ri;
      csound2 = gammag * pl / rl;
      a_star2 = (csound2 + bsquared);
      term = sqrt(a_star2 * a_star2 - 4 * csound2 * calfven2);
      cfast2 = 0.5 * (a_star2 + term);
      cfast = sqrt(cfast2);

      speed1 = fabs(ul) + cfast;
      *maxspeed = (max)(speed1, *maxspeed);

      calfven2 = bv * bv * ri;
      term = sqrt(a_star2 * a_star2 - 4 * csound2 * calfven2);
      cfast2 = 0.5 * (a_star2 + term);
      cfast = sqrt(cfast2);
      speed2 = fabs(vl) + cfast;
      *maxspeed = (max)(speed2, *maxspeed);

      calfven2 = bw * bw * ri;
      term = sqrt(a_star2 * a_star2 - 4 * csound2 * calfven2);
      cfast2 = 0.5 * (a_star2 + term);
      cfast = sqrt(cfast2);
      speed3 = fabs(wl) + cfast;
      *maxspeed = (max)(speed3, *maxspeed);


      /*
         cout << ii
         << " " << jj
         << " " << rl
         << " " << px
         << " " << py
         << " " << pl
         << endl;
         cout << *maxspeed << endl;
         if ( *maxspeed > 700)
         exit(0);
       */
    }
  }

  return 0;

}
