#ifndef __PeningTrap_hpp__
#define __PenningTrap_hpp__

#include <vector>
#include <armadillo>
#include "particle.hpp"

class PenningTrap 
{
public:
  // Constructor
  PenningTrap(double B0, double V0, double d, std::vector<Particle> p);

  double B0;               // Magnetic field strength
  double V0;               // Applied potential
  double d;                // Charismatic dimension
  std::vector<Particle> p; // Particled

  // Methods
  arma::vec external_E_field(arma::vec r); 
  arma::vec external_B_field(arma::vec r); 
  arma::vec force_particle(int i, int j);


};

#endif