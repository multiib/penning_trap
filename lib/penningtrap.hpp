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

  double f;
  double w_V;

  // Methods
  arma::vec external_E_field(arma::vec r); 

  arma::vec external_B_field(arma::vec r); 
  arma::vec force_particle(int i, int j);
  arma::vec total_force_external(int i);

  arma::vec total_force_particles(int i, bool interact);
  arma::vec total_force(int i, bool interact);

  void evolve_RK4(double dt, bool interact);

  void evolve_forward_Euler(double dt, bool interact);

  int count_particles(void);
  void set_amp_freq(double amp, double freq);

  arma::vec external_E_field(arma::vec r, double t);
  arma::vec total_force_external(int i, double t);
  arma::vec total_force(int i, bool interact, double t);
  void evolve_RK4(double dt, bool interact, double t);
};

#endif