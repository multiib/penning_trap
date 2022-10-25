#ifndef __Analytical_hpp__
#define __Analytical_hpp__

#include <armadillo>

arma::mat analytical(
	arma::vec r_in,
	arma::vec v_in,
	double B0,
	double V0,
	double d,
	double m,
	double q,
	double n,
	double t);




#endif