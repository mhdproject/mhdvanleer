  /* Two-dimensional Upwind/Donor-cell code - taken from Hirsch II */
/* Second order Runge Kutta - uses halfstep */
#include "main.h"
#define JET
//#undef JET
#ifndef JET
#define SHOCKTUBE
#endif

#define VERBOSE_OUTPUT
//#define DEBUG_HALFSTEP
#define SECOND_ORDER_TIME
//#undef SECOND_ORDER_TIME
#undef DEBUG
int ne = 0;
int nx = 0;
int ny = 0;
int nz = 0;
double delta_x = 1.0;

double gammag = 1.666666666666667;
double pmin = 0.01;

double max (double a, double b);

int
main (int argc, char **argv)
{
  /* Allocate a 2d array */
  ofstream fout;
  ofstream gout;
  ofstream fin;
  ofstream logfile;
  ifstream init;


  int inject_jet=0;


  int ii = 0, jj = 0, kk = 0;
  int rc;
  int timestep = 0;
  int maxstep = 0;
  int printtime = 10;
  int second = 0;
  double max_speed = 0.0;
  double maxtime = 0.0;
  double *maximumspeed = &max_speed;
  double time = 0.0;
  double delta_t = 0.0;
  double dtodx = 0.0;
  double del = 0.0;
  double delh = 0.0;
  double cfl = 0.80;
  double gammam1 = gammag - 1;
  double fn[8];
  zone maxvar;
  zone minvar;
#ifdef DEBUG
  double px, py, et, ri, rl, ul, vl, ke, pl, al;
#endif /* DEBUG */

  clock_t start, end;
  double elapsed;
  double den_temp[4];
  string probtype;


  start = clock ();


  gammag = 5.0 / 3.0;
//  gammag =1.4;
//  
  /*
     init.open ("input/init.dat");
     init >> nx >> ny;
     init >> maxstep;
     init >> cfl;
     init.close ();
   */
  logfile.open ("output/roe.log");
  logfile.close ();
  logfile.open ("output/falle.log");
  logfile.close ();
  logfile.open ("logfile.txt");
  logfile.close ();
  ne = 8;


  if (argc > 1)
    {
      init.open (argv[1]);
    }
  else
    {
      init.open ("input/gaz1");
    }

  init >> nx;
  init.ignore (256, '\n');
  init >> ny;
  init.ignore (256, '\n');
  init >> maxstep;
  init.ignore (256, '\n');
  init >> cfl;
  init.ignore (256, '\n');
  init >> printtime;
  init.ignore (256, '\n');
  init >> delta_x;
  init.ignore (256, '\n');
  init >> gammag;
  init.ignore (256, '\n');
  init >> maxtime;
  init.ignore (256, '\n');
  init >> probtype;
  init.close ();


  gammam1 = gammag - 1;


  nz = 1;
  Array3D < zone > grid (nx, ny, nz);
  Array3D < zone > gridh (nx, ny, nz);
  Array3D < zone > gridn (nx, ny, nz);
  Array3D < zone > fx (nx, ny, nz);
  Array3D < zone > fy (nx, ny, nz);
  Array3D < zone > xResState (nx, ny, nz);
  Array3D < zone > yResState (nx, ny, nz);

  delta_x = 1.0 / nx;

  cout << "\t\t\t2D Roe Solver Code" << endl;
  cout << "\t\t\tVersion 2.0" << endl;
  cout << "\t\t\tNo of steps =  " << maxstep << endl;
  cout << "\t\t\tNX    =  " << nx << endl;
  cout << "\t\t\tNY    =  " << ny << endl;
  cout << "\t\t\tNE    =  " << ne << endl;
  cout << "\t\t\tCFL   =  " << cfl << endl;
  cout << "\t\t\tdel_x =  " << delta_x << endl;
  cout << "\t\t\tGamma =  " << gammag << endl;
#ifdef SECOND_ORDER_TIME
  cout << "\t\t\t2nd Order Time and Space " << endl;
#endif /* SECOND_ORDER_TIME */
#ifdef LAPIDUS_VISCOSITY
  cout << "\t\t\tLapidus Viscosity " << endl;
#endif /* LAPIDUS_VISCOSITY */
  cout << endl;
  cout << endl;

  // Set up initial values of conserved variables on the grid */


  if (probtype == "Shock")
    {
      if (argc > 1)
	{
	  rc = initialise (argv[1], grid, &maxstep, &cfl);
	}
      else
	{
	  rc = initialise ("input/gaz1", grid, &maxstep, &cfl);
	}

    }

  else if (probtype == "Jet")
    {
      if (argc > 1)
	{
	  rc = initialise_jet (argv[1], grid, &maxstep, &cfl);
	}
      else
	{
	  rc = initialise_jet ("input/input.jet", grid, &maxstep, &cfl);
	}
    }

  else if (probtype == "Blast")
    {
      if (argc > 1)
	{
	  rc = initialise_blast (argv[1], grid, &maxstep, &cfl);
	}
      else
	{
	  rc = initialise_blast ("input/input.jet", grid, &maxstep, &cfl);
	}
    }



/* Print out the initial array */
  fin.open ("infile.txt");
  for (ii = 0; ii < nx; ii++)
    {
      for (jj = 0; jj < ny; jj++)
	{
	  fin << ii
	    << " " << jj
	    << " " << grid[ii][jj][kk] _MASS
	    << " " << grid[ii][jj][kk] _MOMX
	    << " " << grid[ii][jj][kk] _MOMY
	    << " " << grid[ii][jj][kk] _MOMZ
	    << " " << grid[ii][jj][kk] _ENER
	    << " " << grid[ii][jj][kk] _B_X
	    << " " << grid[ii][jj][kk] _B_Y
	    << " " << grid[ii][jj][kk] _B_Z << endl;
	}
    }
  fin.close ();


#ifdef DEBUG1
  for (ii = 0; ii < nx; ii++)
    {
      for (jj = 0; jj < ny; jj++)
	{
	  cout << grid[ii][jj][kk] _MASS << "\t";
	}
      cout << endl;
    }

#endif
/* Using a second order in time Runge-Kutta method, advect the array */

	  rc = output (grid, fx, fy, 0, "out_2d_");
  for (timestep = 1; timestep < maxstep; timestep++)
    {
      for (int k = 0; k < ne; k++)
	{
	  maxvar.array[k] = 0.;
	  minvar.array[k] = 999.;
	}
      /* Set the maximum wave speed to zero */
      *maximumspeed = 0;
      /* Determine the maximum wave speed for use in
       * the time step */
      rc = maxspeed (grid, maximumspeed);

      /* Determine a value for time advance and courant number based on the
       * maximum wave speed */

      del = cfl / *maximumspeed;
      delh = 0.5 * del;
      delta_t = del * delta_x;
		dtodx = 1.0/ *maximumspeed;
      time = time + delta_t;
#ifdef VERBOSE_OUTPUT
      cout << "Tstep= " << setiosflags (ios::scientific) << timestep;
      cout << "\tTime= " << setiosflags (ios::scientific) << time;
      cout << "\tMaxspeed= " << setiosflags (ios::
					     scientific) << *maximumspeed;
      cout << "\t dt = " << setiosflags (ios::scientific) << delta_t;
      cout << "\tCFL= " << setiosflags (ios::scientific) << del;
      cout << endl;
#endif /* VERBOSE_OUTPUT */


//  rc = output ( grid, fx, fy, timestep, "oldg_2d_");

#ifdef SECOND_ORDER_TIME
      jj = 0;


#ifdef TWODIM
      for (jj = 2; jj < ny - 1; jj++)
#endif /* TWODIM */
	{
	  for (ii = 2; ii < nx - 1; ii++)
	    {
	      rc =
		flux (grid, fx[ii][jj][kk].array, xResState[ii][jj][kk].array, dtodx,
		      ii, jj, timestep, 1, 0);
#ifdef TWODIM
	      rc =
		flux (grid, fy[ii][jj][kk].array, yResState[ii][jj][kk].array,dtodx,
		      ii, jj, timestep, 2, 0);
#endif /* TWODIM */
	    }
	}

//#ifdef TWODIM
//          for (jj = 2; jj < ny - 2; jj++)
//#endif /* TWODIM */
//        {
      //      for (ii = 2; ii < nx - 2; ii++)
//             {

      rc =
	update (gridh, grid, fx, fy, xResState, yResState, delh, ii, jj,
		timestep, grid, delta_t, 0);

//             }
//        }

      /* Boundary Conditions */
      rc = boundary (gridh, inject_jet);
      /* End Boundary Conditions */
#ifdef DEBUG_HALFSTEP
      if (timestep % printtime == 0)
	{
	  rc = output (gridh, fx, fy, timestep, "hout_2d_");
	}
#endif /* DEBUG HALFSTEP */


#ifdef TWODIM
      for (jj = 2; jj < ny - 1; jj++)
#endif /* TWODIM */
	{
	  for (ii = 2; ii < nx - 1; ii++)
	    {
	      rc =
		flux (gridh, fx[ii][jj][kk].array,
		      xResState[ii][jj][kk].array,dtodx, ii, jj, timestep, 1, 1);
#ifdef TWODIM
	      rc =
		flux (gridh, fy[ii][jj][kk].array,
		      yResState[ii][jj][kk].array,dtodx, ii, jj, timestep, 2, 1);
#endif /* TWODIM */
	    }
	}

      rc =
	update (gridn, grid, fx, fy, xResState, yResState, del, ii, jj,
		timestep, gridh, delta_t, 1);
#ifdef TWODIM
      for (jj = 2; jj < ny - 2; jj++)
#endif /* TWODIM */
	{
	  for (ii = 2; ii < nx - 2; ii++)
	    {
	      for (int k = 0; k < ne; k++)
		{
		  maxvar.array[k] =
		    (maxvar.array[k] >
		     gridn[ii][jj][kk].array[k] ? maxvar.
		     array[k] : gridn[ii][jj][kk].array[k]);
		  minvar.array[k] =
		    (minvar.array[k] <
		     gridn[ii][jj][kk].array[k] ? minvar.
		     array[k] : gridn[ii][jj][kk].array[k]);
		}

	    }
	}
      grid = gridn.copy ();
//      rc = output ( gridh, fx, fy, timestep, "hout_2d_");
#else /* First order */


#ifdef TWODIM
      for (jj = 2; jj < ny - 1; jj++)
#endif /* TWODIM */
	{
	  for (ii = 2; ii < nx - 1; ii++)
	    {
	      rc =
		flux (grid, fx[ii][jj][kk].array, xResState[ii][jj][kk].array,dtodx,
		      ii, jj, timestep, 1, 0);
#ifdef TWODIM
	      rc =
		flux (grid, fy[ii][jj][kk].array, yResState[ii][jj][kk].array,dtodx,
		      ii, jj, timestep, 2, 0);
#endif /* TWODIM */
	    }
	}

      rc =
	update (gridn, grid, fx, fy, xResState, yResState, del, ii, jj,
		timestep, grid, delta_t, 0);
#ifdef TWODIM
      for (jj = 2; jj < ny - 2; jj++)
#endif /* TWODIM */
	{
	  for (ii = 2; ii < nx - 2; ii++)
	    {

	      for (int k = 0; k < ne; k++)
		{
		  maxvar.array[k] =
		    (maxvar.array[k] >
		     gridn[ii][jj][kk].array[k] ? maxvar.
		     array[k] : gridn[ii][jj][kk].array[k]);
		  minvar.array[k] =
		    (minvar.array[k] <
		     gridn[ii][jj][kk].array[k] ? minvar.
		     array[k] : gridn[ii][jj][kk].array[k]);
		}
	    }
	}
      grid = gridn.copy ();
#endif
      /* Boundary Conditions */
      rc = boundary (grid, inject_jet);
      /* End Boundary Conditions */

      cout << setiosflags (ios::fixed);
      for (int k = 0; k < ne; k++)
	{
	  cout << k +
	    1 << " " << maxvar.array[k] << " " << minvar.array[k] << endl;
	}


      if (timestep % printtime == 0)
	{
	  // cout << "outputting " << endl;
	  rc = output (grid, fx, fy, timestep, "out_2d_");
	}
      if (time > maxtime)
	{
	  rc = output (grid, fx, fy, timestep, "out_2d_");
	  break;
	}

      /* ------------- */
      /* End of loop through timestep */
    }
  /* ------------- */





  end = clock ();
  elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
  cout << "\n\n Elapsed time = " << elapsed << " sec\n\n" << endl;


  return 0;
}

double
max (double a, double b)
{
//      cout << a << " " << b << endl;
  if (a >= b)
    {
      return a;
    }
  else
    {
      return b;
    }
}
