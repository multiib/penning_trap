#include "../lib/penningtrap.hpp"
#include <vector>
#include <armadillo>
#include <cmath>

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
    if (arma::norm(r) > d)
    {
        return arma::vec("0.0 0.0 0.0");
    }
    
    arma::vec E = arma::vec(3);

    E(0) =  (V0/(d*d))   * r(0);
    E(1) =  (V0/(d*d))   * r(1);
    E(2) = -(2*V0/(d*d)) * r(2);

    return E;
}



// External magnetic field
arma::vec PenningTrap::external_B_field(arma::vec r)
{
    if (arma::norm(r) > d)
    {
        return arma::vec("0.0 0.0 0.0");
    }

    arma::vec B = arma::vec(3);

    B(0) = 0;
    B(1) = 0;
    B(2) = B0;

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
arma::vec PenningTrap::total_force_particles(int i, bool interact)
{
    arma::vec total = arma::vec("0.0 0.0 0.0");

    if (interact == false)
    {
        return total;
    }

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
arma::vec PenningTrap::total_force(int i, bool interact)
{

    
    return total_force_external(i) + total_force_particles(i, interact);
}



// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt, bool interact)
{
    // Create a temp copy of the the particles contained in PenningTrap::particles
    std::vector<Particle> p_old = p;

    // Collections of arma::vec's to store the RK4 k-values for each particle
    std::vector<arma::vec> k1_r(p.size());
    std::vector<arma::vec> k1_v(p.size());
    std::vector<arma::vec> k2_r(p.size());
    std::vector<arma::vec> k2_v(p.size());
    std::vector<arma::vec> k3_r(p.size());
    std::vector<arma::vec> k3_v(p.size());
    std::vector<arma::vec> k4_r(p.size());
    std::vector<arma::vec> k4_v(p.size());

    arma::vec acc = arma::vec(3);
    arma::vec vel = arma::vec(3);


    // K1
    for (int i = 0; i <p.size(); i++)
    {
        acc = total_force(i, interact)/p[i].m;
        vel = p[i].v;

        
        k1_r[i] = dt*vel;
        k1_v[i] = dt*acc;

    }

    // K2
    for (int i = 0; i <p.size(); i++)
    {

        
        p[i].r = p[i].r + k1_r[i]/2;
        p[i].v = p[i].v + k1_v[i]/2;

        acc = total_force(i, interact)/p[i].m;
        vel = p[i].v;

        k2_r[i] = dt*vel;
        k2_v[i] = dt*acc;

    }

    // K3
    for (int i = 0; i <p.size(); i++)
    {
        p[i].r = p[i].r + k2_r[i]/2;
        p[i].v = p[i].v + k2_v[i]/2;

        acc = total_force(i, interact)/p[i].m;
        vel = p[i].v;

        k3_r[i] = dt*vel;
        k3_v[i] = dt*acc;

    }

    // K4
    for (int i = 0; i <p.size(); i++)
    {
        p[i].r = p[i].r + k3_r[i];
        p[i].v = p[i].v + k3_v[i];

        acc = total_force(i, interact)/p[i].m;
        vel = p[i].v;

        k4_r[i] = dt*vel;
        k4_v[i] = dt*acc;

    }

    // Final update
    for (int i = 0; i <p.size(); i++)
    {
        p[i].r = p_old[i].r + (k1_r[i] + 2*k2_r[i] + 2*k3_r[i] + k4_r[i])/6;
        p[i].v = p_old[i].v + (k1_v[i] + 2*k2_v[i] + 2*k3_v[i] + k4_v[i])/6;
    }
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt, bool interact)
{
    std::vector<Particle> p_old = p;

    for (int i = 0;i <p.size(); i++)
    {

        arma::vec acc = arma::vec(3);
        
        acc = total_force(i, interact)/p[i].m;
        p[i].v = p_old[i].v + dt * acc;
        p[i].r = p_old[i].r + dt * p[i].v;

    }

}

// Count how many of the particles are still inside the trap region
int PenningTrap::count_particles(void)
{
    int count = 0;
    for (int i = 0; i < p.size(); i++)
    {
        if (arma::norm(p[i].r) < d)
        {
            count++;
        }
    }
    return count;
}

void PenningTrap::set_amp_freq(double amp, double freq)
{
    f = amp;
    w_V = freq;
}

// TIME DEPENDENT OVERLOADED METHODS


// External time dependent electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r, double t)
{
    if (arma::norm(r) > d)
    {
        return arma::vec("0.0 0.0 0.0");
    }


    double V0_t = V0*(1 + f*cos(w_V*t));

    arma::vec E = arma::vec(3);

    E(0) =  (V0_t/(d*d))   * r(0);
    E(1) =  (V0_t/(d*d))   * r(1);
    E(2) = -(2*V0_t/(d*d)) * r(2);

    return E;
}

// The total force on particle_i from the external fields TIME DEPENDENT
arma::vec PenningTrap::total_force_external(int i, double t)
{


    return (p[i].q*external_E_field(p[i].r, t)
           + arma::cross(p[i].v, external_B_field(p[i].r)));
}

// The total force on particle_i from both external fields and other particles TIME DEPENDENT
arma::vec PenningTrap::total_force(int i, bool interact, double t)
{

    return total_force_external(i, t) + total_force_particles(i, interact);
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt, bool interact, double t)
{
    // Create a temp copy of the the particles contained in PenningTrap::particles
    std::vector<Particle> p_old = p;

    // Collections of arma::vec's to store the RK4 k-values for each particle
    std::vector<arma::vec> k1_r(p.size());
    std::vector<arma::vec> k1_v(p.size());
    std::vector<arma::vec> k2_r(p.size());
    std::vector<arma::vec> k2_v(p.size());
    std::vector<arma::vec> k3_r(p.size());
    std::vector<arma::vec> k3_v(p.size());
    std::vector<arma::vec> k4_r(p.size());
    std::vector<arma::vec> k4_v(p.size());

    arma::vec acc = arma::vec(3);
    arma::vec vel = arma::vec(3);


    // K1
    for (int i = 0; i <p.size(); i++)
    {
        acc = total_force(i, interact, t)/p[i].m;
        vel = p[i].v;

        
        k1_r[i] = dt*vel;
        k1_v[i] = dt*acc;

    }

    // K2
    for (int i = 0; i <p.size(); i++)
    {

        
        p[i].r = p[i].r + k1_r[i]/2;
        p[i].v = p[i].v + k1_v[i]/2;

        acc = total_force(i, interact, t+.5*dt)/p[i].m;
        vel = p[i].v;

        k2_r[i] = dt*vel;
        k2_v[i] = dt*acc;

    }

    // K3
    for (int i = 0; i <p.size(); i++)
    {
        p[i].r = p[i].r + k2_r[i]/2;
        p[i].v = p[i].v + k2_v[i]/2;

        acc = total_force(i, interact, t+.5*dt)/p[i].m;
        vel = p[i].v;

        k3_r[i] = dt*vel;
        k3_v[i] = dt*acc;

    }

    // K4
    for (int i = 0; i <p.size(); i++)
    {
        p[i].r = p[i].r + k3_r[i];
        p[i].v = p[i].v + k3_v[i];

        acc = total_force(i, interact, t+dt)/p[i].m;
        vel = p[i].v;

        k4_r[i] = dt*vel;
        k4_v[i] = dt*acc;

    }

    // Final update
    for (int i = 0; i <p.size(); i++)
    {
        p[i].r = p_old[i].r + (k1_r[i] + 2*k2_r[i] + 2*k3_r[i] + k4_r[i])/6;
        p[i].v = p_old[i].v + (k1_v[i] + 2*k2_v[i] + 2*k3_v[i] + k4_v[i])/6;
    }
}