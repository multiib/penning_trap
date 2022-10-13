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

  double B0;    // Magnetic field strength
  double V0;    // Applied potential
  double d; // Charismatic dimension
  std::vector<Particle> p; // Particled

};

#endif