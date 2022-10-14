#include "lib/penningtrap.hpp"
#include <vector>
#include <armadillo>

// Definition of the constructor
PenningTrap::PenningTrap(double B0, double V0, double d, std::vector<Particle> p)
{
    // Use the input variables (c0, c1) to assign values to the class memeber variables (c0_, c1_)
    B0 = B0;
    V0 = V0;
    d  = d;
    p  = p;

    // Coulomb constant
    double k_e = 1.38935333e5;
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r)
{
    arma::vec E = arma::vec(3);

    E(0) =  (V0/(d*d))   * r(0);
    E(1) =  (V0/(d*d))   * r(1);
    E(2) = -(2*V0/(d*d)) * r(2);

    return E;
}

// External magnetic field
arma::vec PenningTrap::external_B_field(arma::vec r)
{
    arma::vec B = arma::vec(3, 0);

//  B(0) = 0
//  B(1) = 0
    B(2) = B0 * r(2);

    return B;
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j)
{
    arma::vec ri = p[i].r;
    arma::vec rj = p[j].r;
    
    double xi, xj, yi, yj, zi, zj;

    xi = ri(0);
    xj = rj(0);

    yi = ri(1);
    yj = rj(1);

    zi = ri(2);
    zj = rj(2);

    // Normal vector
    double nrm = arma::norm(ri - rj);

}