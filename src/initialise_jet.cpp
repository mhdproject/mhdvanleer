#include "initialise.h"
int
initialise_jet (char *filename, Array3D < zone > grid, int *maxstep,
		double *cfl)
{
  int ii = 0;
  int jj = 0;
  int kk = 0;
  double rr, vx, vy, vz, p, ke, b1, b2, b3;
  double rl, ul, vl, wl, pl, bul, bvl, bwl, al;
  double rho, ur, vr, wr, pr, bur, bvr, bwr, ar;
  ifstream input;


  int tempint;
  double tempdouble;


  double gammam1 = gammag - 1;
  double gammam1i = 1.0 / gammam1;
  double dist;

  double bsquared = 0;
  double sqr4pie = 0;
  double sqr4piei = 0;
  double unused = 0;


  double pie = acos (-1.0);
  sqr4pie = 2.0 * sqrt (pie);
  sqr4piei = 1.0 / sqr4pie;
//  sqr4piei = 1.0 ;


  rl = 1.0;
  ul = 0.0;
  vl = 0.0;
  pl = 6.0;

  rr = 0.1;
  ur = 0.0;
  vr = 0.0;
  pr = 0.6;

  input.open (filename);

  input >> tempint;
  input.ignore (256, '\n');
  input >> tempint;
  input.ignore (256, '\n');
  input >> tempint;
  input.ignore (256, '\n');
  input >> tempdouble;
  input.ignore (256, '\n');
  input >> tempint;
  input.ignore (256, '\n');
  input >> tempdouble;
  input.ignore (256, '\n');
  input >> tempdouble;
  input.ignore (256, '\n');
  input >> tempint;
  input.ignore (256, '\n');
  input.ignore (256, '\n');

  input >> rl;
  input.ignore (256, '\n');
  input >> ul;
  input.ignore (256, '\n');
  input >> vl;
  input.ignore (256, '\n');
  input >> wl;
  input.ignore (256, '\n');
  input >> pl;
  input.ignore (256, '\n');
  input >> bul;
  input.ignore (256, '\n');
  input >> bvl;
  input.ignore (256, '\n');
  input >> bwl;
  input.ignore (256, '\n');

  input >> rr;
  input.ignore (256, '\n');
  input >> ur;
  input.ignore (256, '\n');
  input >> vr;
  input.ignore (256, '\n');
  input >> wr;
  input.ignore (256, '\n');
  input >> pr;
  input.ignore (256, '\n');
  input >> bur;
  input.ignore (256, '\n');
  input >> bvr;
  input.ignore (256, '\n');
  input >> bwr;
  input.ignore (256, '\n');

  input.close ();



  al = sqrt (gammag * pl / rl);
  ar = sqrt (gammag * pr / rr);




  cout << "Ambient state rho= " << rl << "\t Jet State " << rr << endl;
  cout << "Ambient state vx = " << ul << "\t Jet State " << ur << endl;
  cout << "Ambient state vy = " << vl << "\t Jet State " << vr << endl;
  cout << "Ambient state vz = " << wl << "\t Jet State " << wr << endl;
  cout << "Ambient state p  = " << pl << "\t Jet State " << pr << endl;
  cout << "Ambient state Bx = " << bul << "\t Jet State " << bur << endl;
  cout << "Ambient state By = " << bvl << "\t Jet State " << bvr << endl;
  cout << "Ambient state Bz = " << bwl << "\t Jet State " << bwr << endl;
  cout << "Ambient state A  = " << al << "\t Jet State " << ar << endl;
  cout << "Jet Mach number  = " << vr / ar << endl;
  cout << endl;


  /* Initialise the grid array with a shock tube problem  */
  for (ii = 0; ii < nx; ii++)
    {
      for (jj = 0; jj < ny; jj++)
	{
	  vx = ul;
	  vy = vl;
	  vz = wl;
	  p = pl;
	  rho = rl;
	  b1 = bul * sqr4piei;
	  b2 = bvl * sqr4piei;
	  b3 = bwl * sqr4piei;

//   JET
	  double xdist = ii - nx / 2;
	  double ydist = jj - ny / 2;
	  dist = sqrt ((double) (xdist * xdist + ydist * ydist));
//	  dist = sqrt ((double) (xdist * xdist ));
	  if (dist < 0.1 * nx )
	    // if ( ii +jj < 0.5*nx ) 

	    {
	      vx = ur;
	      vy = vr;
	      vz = wr;
	      p = pr;
	      rho = rr;
	      b1 = bur * sqr4piei;
	      b2 = bvr * sqr4piei;
	      b3 = bwr * sqr4piei;
	    }

	  bsquared = (b1 * b1 + b2 * b2 + b3 * b3);
	  grid[ii][jj][kk] _MASS = rho;
	  grid[ii][jj][kk] _MOMX = rho * vx;
	  grid[ii][jj][kk] _MOMY = rho * vy;
	  grid[ii][jj][kk] _MOMZ = rho * vz;
	  ke = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
	  grid[ii][jj][kk] _ENER = ke + p * gammam1i + 0.5 * bsquared;
	  grid[ii][jj][kk] _B_X = b1;
	  grid[ii][jj][kk] _B_Y = b2;
	  grid[ii][jj][kk] _B_Z = b3;
	  if (grid[ii][jj][kk] _ENER < 0.0)
	    {
	      cout << "Wtf?" << endl;
	      exit (0);
	    }
	}

    }

  return 0;
}
int
initialise_blast (char *filename, Array3D < zone > grid, int *maxstep,
		double *cfl)
{
  int ii = 0;
  int jj = 0;
  int kk = 0;
  double rr, vx, vy, vz, p, ke, b1, b2, b3;
  double rl, ul, vl, wl, pl, bul, bvl, bwl, al;
  double rho, ur, vr, wr, pr, bur, bvr, bwr, ar;
  ifstream input;


  int tempint;
  double tempdouble;


  double gammam1 = gammag - 1;
  double gammam1i = 1.0 / gammam1;
  double dist;

  double bsquared = 0;
  double sqr4pie = 0;
  double sqr4piei = 0;
  double unused = 0;


  double pie = acos (-1.0);
  sqr4pie = 2.0 * sqrt (pie);
  sqr4piei = 1.0 / sqr4pie;
//  sqr4piei = 1.0 ;


  rl = 1.0;
  ul = 0.0;
  vl = 0.0;
  pl = 6.0;

  rr = 0.1;
  ur = 0.0;
  vr = 0.0;
  pr = 0.6;

  input.open (filename);

  input >> tempint;
  input.ignore (256, '\n');
  input >> tempint;
  input.ignore (256, '\n');
  input >> tempint;
  input.ignore (256, '\n');
  input >> tempdouble;
  input.ignore (256, '\n');
  input >> tempint;
  input.ignore (256, '\n');
  input >> tempdouble;
  input.ignore (256, '\n');
  input >> tempdouble;
  input.ignore (256, '\n');
  input >> tempint;
  input.ignore (256, '\n');
  input.ignore (256, '\n');

  input >> rl;
  input.ignore (256, '\n');
  input >> ul;
  input.ignore (256, '\n');
  input >> vl;
  input.ignore (256, '\n');
  input >> wl;
  input.ignore (256, '\n');
  input >> pl;
  input.ignore (256, '\n');
  input >> bul;
  input.ignore (256, '\n');
  input >> bvl;
  input.ignore (256, '\n');
  input >> bwl;
  input.ignore (256, '\n');

  input >> rr;
  input.ignore (256, '\n');
  input >> ur;
  input.ignore (256, '\n');
  input >> vr;
  input.ignore (256, '\n');
  input >> wr;
  input.ignore (256, '\n');
  input >> pr;
  input.ignore (256, '\n');
  input >> bur;
  input.ignore (256, '\n');
  input >> bvr;
  input.ignore (256, '\n');
  input >> bwr;
  input.ignore (256, '\n');

  input.close ();



  al = sqrt (gammag * pl / rl);
  ar = sqrt (gammag * pr / rr);




  cout << "Ambient state rho= " << rl << "\t Jet State " << rr << endl;
  cout << "Ambient state vx = " << ul << "\t Jet State " << ur << endl;
  cout << "Ambient state vy = " << vl << "\t Jet State " << vr << endl;
  cout << "Ambient state vz = " << wl << "\t Jet State " << wr << endl;
  cout << "Ambient state p  = " << pl << "\t Jet State " << pr << endl;
  cout << "Ambient state Bx = " << bul << "\t Jet State " << bur << endl;
  cout << "Ambient state By = " << bvl << "\t Jet State " << bvr << endl;
  cout << "Ambient state Bz = " << bwl << "\t Jet State " << bwr << endl;
  cout << "Ambient state A  = " << al << "\t Jet State " << ar << endl;
  cout << "Jet Mach number  = " << vr / ar << endl;
  cout << endl;


  /* Initialise the grid array with a shock tube problem  */
  for (ii = 0; ii < nx; ii++)
    {
      for (jj = 0; jj < ny; jj++)
	{
	  vx = ul;
	  vy = vl;
	  vz = wl;
	  p = pl;
	  rho = rl;
	  b1 = bul * sqr4piei;
	  b2 = bvl * sqr4piei;
	  b3 = bwl * sqr4piei;

//   JET
	  double xdist = ii - nx / 2;
	  double ydist = jj - ny / 2;
	  dist = sqrt ((double) (xdist * xdist + ydist * ydist));
//	  dist = sqrt ((double) (xdist * xdist ));
	  if (dist < 0.1 * nx )
	    // if ( ii +jj < 0.5*nx ) 

	    {
	      vx = ur;
	      vy = vr;
	      vz = wr;
	      p = pr;
	      rho = rr;
	      b1 = bur * sqr4piei;
	      b2 = bvr * sqr4piei;
	      b3 = bwr * sqr4piei;
	    }

	  bsquared = (b1 * b1 + b2 * b2 + b3 * b3);
	  grid[ii][jj][kk] _MASS = rho;
	  grid[ii][jj][kk] _MOMX = rho * vx;
	  grid[ii][jj][kk] _MOMY = rho * vy;
	  grid[ii][jj][kk] _MOMZ = rho * vz;
	  ke = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
	  grid[ii][jj][kk] _ENER = ke + p * gammam1i + 0.5 * bsquared;
	  grid[ii][jj][kk] _B_X = b1;
	  grid[ii][jj][kk] _B_Y = b2;
	  grid[ii][jj][kk] _B_Z = b3;
	  if (grid[ii][jj][kk] _ENER < 0.0)
	    {
	      cout << "Wtf?" << endl;
	      exit (0);
	    }
	}

    }

  return 0;
}
