#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

class Particle 
{
public:
  // Constructor
  Particle(double charge, double mass, arma::vec position, arma::vec velocity);

  double q;    // Charge of particle
  double m;    // Mass off particle
  arma::vec r; // Initial position
  arma::vec v; // Initial velocity

};

#endif