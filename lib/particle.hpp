#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

class Particle 
{
public:
  // Constructor
  Particle(double q, double m, arma::vec r, arma::vec v);

  double q;    // Charge
  double m;    // Mass
  arma::vec r; // Position
  arma::vec v; // Velocity

};

#endif