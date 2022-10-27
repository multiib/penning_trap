#include "../lib/particle.hpp"

// Definition of the constructor
Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity)
{
  q = charge;   // Charge of particle
  m = mass;     // Mass of particle
  r = position; // Initial position
  v = velocity; // Initial velocity
}