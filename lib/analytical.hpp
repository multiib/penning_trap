#ifndef __Analytical_hpp__
#define __Analytical_hpp__

#include <armadillo>

arma::mat analytical(
	arma::vec r_in, // Initial position (y value has to be 0)
	arma::vec v_in, // Initial velocity (x and z value has to be 0)
	double B0,      // Magnetic field strength
	double V0,      // Potential applied to the electrodes
	double d,       // Characteristic dimension
	double m,       // Mass of particle
	double q,       // Charge of particle
	double n,       // Number of steps
	double t);      // Time of simulation




#endif