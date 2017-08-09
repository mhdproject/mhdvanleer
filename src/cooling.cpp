#include "cooling.h"
#define MOLCOOL
//#undef MOLCOOL


int
RadiativeCooling::cooling(const double *zone, double *Lcooling, double dt) {



// Establish a system of units to work out the cooling
//     Steve's system of units
//
//     LU=10^15cm    ! YSO jet radius
//     VU=10^7 cm/s  ! 100km/s : YSO jet velocity
//     DU=10^-22     ! 120 cm^-3 100% ionised atomic H: typical ambient density
//
//     temperatureU=10^8s      ! 3.17 years
//     MU=10^23g     ! 5*10-11 solar masses
//     EU=10^37ergs  ! 26*10^3 solar luminosity integrated over 1 second
//     PU=10-8barye  ! 60 times pressuresure of reference ISM at 10^4K
//
//     proton mass= 1.67*10-47 MU
//     proton mass only used in deriving number density
//     number density only used in: 1) temperature=p/n*k 2) Lnet=n^2*LN
//     absorb 10^+47 into k and LN
//


  rho = zone[0];
  px = zone[1];
  py = zone[2];
  et = zone[3];

  //     cout << "cool:" ;

  // Work out temperature from zone value
  //
  rhoi = 1.0 / rho;
  velx = px * rhoi;
  vely = py * rhoi;
  velx2 = velx * velx;
  vely2 = vely * vely;
  ke = 0.5 * rho * (velx2 + vely2);
  pressure = et - ke;
  pressure = pressure * (gammam1);

  nt = 2 * mpi * rho;
  nt2 = nt * nt;
  temperature = ki * pressure / nt;
  real_temp = ki * pressure / nt;
  temperature = log10(temperature);
  // set up initial values 

// read cooling table


// Set the minimum subcycling step
  /* initialisation of cooling variables */
  lowest_temperature = atomic_temp_tab[0];    // lowest tabulated temp (LOG)
  subttot = 0.0;
  firststep = 1;

// set the initial energy
  e_init = et - ke;
  et = et - ke;

#ifdef MOLCOOL
  if (temperature < 3000) {
    chi = 0;
  } else if (temperature > 3000 && temperature < 12600) {
    chi = (min)((temperature - 3000) / 7000, 0.9);
  } else {
    chi = 0;
  }
  nH2 = nt * (1 - chi);
  nh = nt * chi;
  rc = molcool(real_temp, nH2, nh, &molrate);
  assert(rc == 0);
#endif /* MOLCOOL */
  if (temperature <= lowest_temperature) {
    rate = 0.0;
  } else {
//    need to use lookup table here
    tabfind(temperature, &rate, atomic_temp_tab);    // absorbed 10+94
  }
#ifdef MOLCOOL
  rate = rate + molrate;
#endif /* MOLCOOL */


//    initialising substep size and initial energy loss rate for all atoms
  eloss = nt2 * rate;
// loop around the substeps 
  *Lcooling = 0;

  counter = 0;
  // Substepping
  while (subttot < dt && temperature > lowest_temperature) {
    double w,  elosstot;
    w = dt - subttot;
    y = 0.05 * et / abs(eloss);    //    max timestep for 5% change in e
    subdt = (min)(w, y);
    elosstot = eloss * subdt;
//        write(*,*) counter,temperature,'precool:t', e/elosstot
    et = et - elosstot;

#ifdef MOLCOOL
    if (temperature < 3000) {
      chi = 0;
    } else if (temperature > 3000 && temperature < 12600) {
      chi = (min)((temperature - 3000) / 7000, 0.9);
    } else {
      chi = 0;
    }
    nH2 = nt * (1 - chi);
    nh = nt * chi;
    rc = molcool(real_temp, nH2, nh, &molrate);
    assert(rc == 0);
//   rc = molcool (pow(10,temperature),nH2, nh  ,&rate_temp );
#endif /* MOLCOOL */
    if (temperature <= lowest_temperature) {
      eloss = 0.0;
    } else {
      tabfind(temperature, &rate, atomic_temp_tab);
#ifdef MOLCOOL
      rate = rate + molrate;
#endif /* MOLCOOL */
      eloss = nt2 * rate;
    }

    counter = counter + 1;
    subttot = subttot + subdt;

    if (firststep == 0) {
      firststep = 0;
    }

  }                //     end of substepping

//

  if (counter > 0) {
    de = et - e_init;        //     energy correction
    *Lcooling = de;
  }

  return 0;
}
