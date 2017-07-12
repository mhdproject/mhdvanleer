/* $Id: godunov_mhd.h,v 1.2 2005/07/25 13:44:24 gmurphy Exp $ */
//#include <math.h>
//#include <stdio.h>
//#include "burger.h"
#include "global.h"
#include "falle.h"
int      godunov_mhd (zone ** grid, zone ** flux1, zone ** flux2,
		  int height, int width, int jj, double  * time);
