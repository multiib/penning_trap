#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <armadillo>
#include <cmath>

#include "../lib/particle.hpp"
#include "../lib/penningtrap.hpp"
#include "../lib/particle.hpp"
#include "../lib/analytical.hpp"

int simulate(
    PenningTrap trap,
    double time,         
    double steps,
    std::string outfile, // Name of outfile
    bool interact,       // Columb interaction switch
    bool FE,
    bool timedependent);            // Forward Euler switch

std::vector<Particle> fill_rand(int n_particles, double d, double q, double m);


// In the main function we find all the simulations and experiments ran in this
// project. The code is written to easily do own experiments.


int main()
{
    // Seed for reproducibility
    arma::arma_rng::set_seed(13141);

    // Simulation configurations and trap properties
    double time = 50.0; // µs
    int steps = 1e6;
    double B0 = 9.65e1; // 1.00 T = u/((µs)**2*e)
    double V0 = 2.41e6; // 25.o mV = u(µm)**2/((µs)**2*e)
    double d = 500; // µm
    bool interact;
    std::string outfile;

    // Ca+ properties
    double q = 1.0;    // Charge +1 e
    double m = 40.078; // Mass of Ca+


    // MAIN SIMULATIONS

    // Particle 1
    arma::vec r1 = arma::vec("20.0 0.0 20.0"); // µm
    arma::vec v1 = arma::vec("0.0 25.0 0.0"); // µm/µs
    Particle p1 = Particle(q,m,r1,v1);

    // Particle 2
    arma::vec r2 = arma::vec("25.0 25.0 0.0"); // µm
    arma::vec v2 = arma::vec("0.0 40.0 5.0"); // µm/µs
    Particle p2 = Particle(q,m,r2,v2);

    // Vector of particles, double simulation
    std::vector<Particle> p_double;
    p_double.push_back(p1);
    p_double.push_back(p2);

    // Double particle simulation
    PenningTrap t_double = PenningTrap(B0, V0, d, p_double);
    interact = false;
    outfile = "./res/out/double.bin";
    simulate(t_double, time, steps, outfile, interact, false, false);

    // Double particle simulation with interactions
    interact = true;
    outfile = "./res/out/double_interactions.bin";
    simulate(t_double, time, steps, outfile, interact, false, false);



    // ERROR SIMULATIONS

    std::vector<Particle> p_single;
    p_single.push_back(p1);

    PenningTrap t_single = PenningTrap(B0, V0, d, p_single);

    // n = 4000
    arma::mat a1 = analytical(r1,v1,B0,V0,d,m,q,4000.0,time);
    outfile = "./res/out/re_AN/4000.bin";
    a1.save(outfile, arma::arma_binary);
    outfile = "./res/out/re_RK4/4000.bin";
    simulate(t_single, time, 4000, outfile, false, false, false);
    outfile = "./res/out/re_FE/4000.bin";
    simulate(t_single, time, 4000, outfile, false, true, false);

    // n = 8000
    arma::mat a2 = analytical(r1,v1,B0,V0,d,m,q,8000.0,time);
    outfile = "./res/out/re_AN/8000.bin";
    a2.save(outfile, arma::arma_binary);
    outfile = "./res/out/re_RK4/8000.bin";
    simulate(t_single, time, 8000, outfile, false, false, false);
    outfile = "./res/out/re_FE/8000.bin";
    simulate(t_single, time, 8000, outfile, false, true, false);

    // n = 16000
    arma::mat a3 = analytical(r1,v1,B0,V0,d,m,q,16000.0,time);
    outfile = "./res/out/re_AN/16000.bin";
    a3.save(outfile, arma::arma_binary);
    outfile = "./res/out/re_RK4/16000.bin";
    simulate(t_single, time, 16000, outfile, false, false, false);
    outfile = "./res/out/re_FE/16000.bin";
    simulate(t_single, time, 16000, outfile, false, true, false);

    // n = 1600
    arma::mat a4 = analytical(r1,v1,B0,V0,d,m,q,32000.0,time);
    outfile = "./res/out/re_AN/32000.bin";
    a4.save(outfile, arma::arma_binary);
    outfile = "./res/out/re_RK4/32000.bin";
    simulate(t_single, time, 32000, outfile, false, false, false);
    outfile = "./res/out/re_FE/32000.bin";
    simulate(t_single, time, 32000, outfile, false, true, false);


    // RESONANCE SIMULATIONS

    double w_V_start, w_V_end, step_size;
    int N, n_particles;
    arma::vec w_V_vec;
    arma::mat M_out;

    // Amplitudes
    arma::vec f_vec = arma::vec("0.1 0.4 0.7");


    // Main resonance scan

    // Frequencies
    w_V_start = 0.2;  // MHz
    w_V_end = 2.5;    // MHz
    step_size = 0.02;

    N = std::round((w_V_end - w_V_start)/step_size + 1);

    n_particles = 100;


    w_V_vec = arma::linspace(w_V_start,w_V_end,N);
    M_out = arma::mat(3,w_V_vec.size());
    for (int i = 0; i<3; i++)
    {
        for (int j = 0; j<w_V_vec.size(); j++)
        {
            std::vector<Particle> particles = fill_rand(n_particles,d,q,m);
            PenningTrap t_rand = PenningTrap(B0, V0, d, particles);
            t_rand.set_amp_freq(f_vec(i), w_V_vec(j));
            std::cout<<j<<std::endl;
            M_out(i,j) = simulate(t_rand, 500, 1e4, "NO_OUTPUT", false, false, true);
        }
    }
    outfile = "./res/out/rand.bin";
    M_out.save(outfile, arma::arma_binary);


    // Zoom no interaction

    // Frequencies
    w_V_start = 1.1;  // MHz
    w_V_end = 1.7;    // MHz
    step_size = 0.005;

    N = std::round((w_V_end - w_V_start)/step_size + 1);

    n_particles = 100;


    w_V_vec = arma::linspace(w_V_start,w_V_end,N);
    M_out = arma::mat(3,w_V_vec.size());
    for (int i = 0; i<3; i++)
    {
        for (int j = 0; j<w_V_vec.size(); j++)
        {
            std::vector<Particle> particles = fill_rand(n_particles,d,q,m);
            PenningTrap t_rand = PenningTrap(B0, V0, d, particles);
            t_rand.set_amp_freq(f_vec(i), w_V_vec(j));
            std::cout<<j<<std::endl;
            M_out(i,j) = simulate(t_rand, 500, 1e4, "NO_OUTPUT", true, false, true);
        }
    }
    outfile = "./res/out/rand_Z_F.bin";
    M_out.save(outfile, arma::arma_binary);


    // Zoom with interaction

    // Frequencies
    w_V_start = 1.1;  // MHz
    w_V_end = 1.7;    // MHz
    step_size = 0.005;

    N = std::round((w_V_end - w_V_start)/step_size + 1);

    n_particles = 100;

    w_V_vec = arma::linspace(w_V_start,w_V_end,N);
    M_out = arma::mat(3,w_V_vec.size());
    for (int i = 0; i<3; i++)
    {
        for (int j = 0; j<w_V_vec.size(); j++)
        {
            std::vector<Particle> particles = fill_rand(n_particles,d,q,m);
            PenningTrap t_rand = PenningTrap(B0, V0, d, particles);
            t_rand.set_amp_freq(f_vec(i), w_V_vec(j));
            std::cout<<j<<std::endl;
            M_out(i,j) = simulate(t_rand, 500, 1e4, "NO_OUTPUT", true, false, true);
        }
    }
    outfile = "./res/out/rand_Z_T.bin";
    M_out.save(outfile, arma::arma_binary);

    return 0;
}


int simulate(
    PenningTrap trap,
    double time,
    double steps,
    std::string outfile,
    bool interact,
    bool FE,
    bool timedependent)
{
    // Timestep
    double dt = time/steps;

    // Time vec
    arma::mat t_out = arma::mat(steps, trap.p.size());

    // Position vec
    arma::mat rx_out = arma::mat(steps, trap.p.size());
    arma::mat ry_out = arma::mat(steps, trap.p.size()); 
    arma::mat rz_out = arma::mat(steps, trap.p.size());

    // Velocity vec
    arma::mat vx_out = arma::mat(steps, trap.p.size()); 
    arma::mat vy_out = arma::mat(steps, trap.p.size()); 
    arma::mat vz_out = arma::mat(steps, trap.p.size()); 


    // Out matrix (structure of matrix is expained in the README file)
    arma::mat M_out = arma::mat(steps, 7*trap.p.size());

    // Initial values
    for (int j = 0; j<trap.p.size(); j++)
    {
        M_out(0, j*7+0) = 0.0;
        M_out(0, j*7+1) = trap.p[j].r(0);
        M_out(0, j*7+2) = trap.p[j].r(1);
        M_out(0, j*7+3) = trap.p[j].r(2);
        M_out(0, j*7+4) = trap.p[j].v(0);
        M_out(0, j*7+5) = trap.p[j].v(1);
        M_out(0, j*7+6) = trap.p[j].v(2);
    }

    if (timedependent)
    {
        // Running simulation
        for (int i = 1; i < steps; i++)
        {
            if (FE)
            {
                trap.evolve_forward_Euler(dt, interact, dt*i);  
            }
            else
            {
                trap.evolve_RK4(dt, interact, dt*i);
            }
            
            for (int j = 0; j<trap.p.size(); j++)
            {
                M_out(i, j*7+0) = i*dt;
                M_out(i, j*7+1) = trap.p[j].r(0);
                M_out(i, j*7+2) = trap.p[j].r(1);
                M_out(i, j*7+3) = trap.p[j].r(2);
                M_out(i, j*7+4) = trap.p[j].v(0);
                M_out(i, j*7+5) = trap.p[j].v(1);
                M_out(i, j*7+6) = trap.p[j].v(2);
            }
        }
    }
    else
    {
        // Running simulation
        for (int i = 1; i < steps; i++)
        {
            if (FE)
            {
                trap.evolve_forward_Euler(dt, interact);  
            }
            else
            {
                trap.evolve_RK4(dt, interact);
            }
            
            for (int j = 0; j<trap.p.size(); j++)
            {
                M_out(i, j*7+0) = i*dt;
                M_out(i, j*7+1) = trap.p[j].r(0);
                M_out(i, j*7+2) = trap.p[j].r(1);
                M_out(i, j*7+3) = trap.p[j].r(2);
                M_out(i, j*7+4) = trap.p[j].v(0);
                M_out(i, j*7+5) = trap.p[j].v(1);
                M_out(i, j*7+6) = trap.p[j].v(2);
            }
        }
    }

    // Save file
    if (outfile != "NO_OUTPUT")
    {
        M_out.save(outfile, arma::arma_binary);
    }

    return trap.count_particles();
}


// Fill particle vector with random particles
std::vector<Particle> fill_rand(int n_particles, double d, double q, double m)
{
    std::vector<Particle> particles;

    for (int i = 0; i<n_particles; i++)
    {
        arma::vec r = arma::vec(3).randn() * 0.1 * d;  // random initial position
        arma::vec v = arma::vec(3).randn() * 0.1 * d;  // random initial velocity
        Particle p_i = Particle(q,m,r,v);

        particles.push_back(p_i);
    }
    return particles;
}

