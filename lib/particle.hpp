#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

class Particle 
{
public:
  // Constructor
  Particle(double charge, double mass, arma::vec position, arma::vec velocity);


  double q;    // Charge
  double m;    // Mass
  arma::vec r; // Position
  arma::vec v; // Velocity

};

#endif