#include "global.h"
int      riemann (
		double * leftstate, 
		double * rightstate, 
		double * roeflux, 
		double * res_state, 
		int time_step, 
		double  * max_speed,
		int idir);
