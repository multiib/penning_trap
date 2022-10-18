#include "../lib/particle.hpp"

// Definition of the constructor
Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity)
{
  // Use the input variables (c0, c1) to assign values to the class memeber variables (c0_, c1_)
  q = charge;
  m = mass;
  r = position;
  v = velocity;



}