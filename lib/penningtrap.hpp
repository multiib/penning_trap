#ifndef __PenningTrap_hpp__
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
  std::vector<Particle> p; // Particles

  // Methods
  arma::vec external_E_field(arma::vec r); 
  arma::vec external_B_field(arma::vec r); 
  arma::vec force_particle(int i, int j);
  arma::vec total_force_external(int i);
  arma::vec total_force_particles(int i);
  arma::vec total_force(int i);
  void evolve_RK4(double dt);
  void evolve_forward_Euler(double dt);

};

#endif