#include "../lib/penningtrap.hpp"
#include <vector>
#include <armadillo>

// Definition of the constructor
PenningTrap::PenningTrap(double mag, double pot, double dist, std::vector<Particle> ps)
{
    // Use the input variables (c0, c1) to assign values
    // to the class memeber variables (c0_, c1_)
    B0 = mag;
    V0 = pot;
    d  = dist;
    p  = ps;
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
    arma::vec B = arma::vec(3);

//  B(0) = 0
//  B(1) = 0
    B(2) = B0 * r(2);

    return B;
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j)
{

    Particle p_i = p[i];
    Particle p_j = p[j];

    arma::vec ri = p_i.r;
    arma::vec rj = p_j.r;
    
    double xi, xj, yi, yj, zi, zj;

    xi = ri(0);
    xj = rj(0);

    yi = ri(1);
    yj = rj(1);

    zi = ri(2);
    zj = rj(2);

    // Coulomb constant
    double k_e = 1.38935333e5;

    // Normal vector
    double nrm = arma::norm(ri - rj);


    // Force vector
    arma::vec F = arma::vec(3);
    double charge_vars = k_e*p_i.q*p_j.q/p_i.m;
    F(0) = charge_vars * ((xi - xj)/(nrm*nrm*nrm));
    F(1) = charge_vars * ((yi - yj)/(nrm*nrm*nrm));
    F(2) = charge_vars * ((zi - zj)/(nrm*nrm*nrm));


    return F;
}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i)
{
    arma::vec total = arma::vec(3);

    for (int j = 0; j<p.size(); j++)
    {
        if (i != j)
        {
            total += force_particle(i, j);
        }
    }
    return total;
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i)
{


    return (p[i].q*external_E_field(p[i].r)
           + arma::cross(p[i].v, external_B_field(p[i].r)));
}



// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i)
{

    return total_force_external(i) + total_force_particles(i);
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
{

}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{

      //y_ip1 = y_i + dt * f_i
    for (int i = 0;i <p.size(); i++)
    {

        arma::vec acc = arma::vec(3);
        
        acc = total_force(i)/p[i].m;
        p[i].v = p[i].v + dt * acc;
        p[i].r = p[i].r + dt * p[i].v;
    }

}