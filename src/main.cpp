#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <armadillo>

#include "../lib/particle.hpp"
#include "../lib/penningtrap.hpp"
#include "../lib/particle.hpp"
#include "../lib/analytical.hpp"



void simulate(
    PenningTrap trap,
    double time,
    double steps,
    std::string outfile,
    bool interact,
    bool FE);

std::vector<Particle> fill_rand(int n_particles, double d, double q, double m);


int main()
{

    arma::arma_rng::set_seed(13141);

    // Simulation configurations
    double time = 50.0; // µs
    int steps = 1e6;



    // Penning Trap properties
    double B0 = 9.65e1; // 1.00 T = u/((µs)**2*e)
    double V0 = 2.41e6; // 25.o mV = u(µm)**2/((µs)**2*e)
    double d = 500; // µm


    // Ca+ properties
    double q = 1.0;    // Charge +1 e
    double m = 40.078; // Mass of Ca+


    // Particle 1
    arma::vec r1 = arma::vec("20.0 0.0 20.0"); // µm
    arma::vec v1 = arma::vec("0.0 25.0 0.0"); // µm/µs
    Particle p1 = Particle(q,m,r1,v1);


    // Particle 2
    arma::vec r2 = arma::vec("25.0 25.0 0.0"); // µm
    arma::vec v2 = arma::vec("0.0 40.0 5.0"); // µm/µs
    Particle p2 = Particle(q,m,r2,v2);


    // Vector of particles, single simulation
    std::vector<Particle> p_single;
    p_single.push_back(p1);


    // Vector of particles, double simulation
    std::vector<Particle> p_double;
    p_double.push_back(p1);
    p_double.push_back(p2);


    // Penning Trap 
    PenningTrap t_single = PenningTrap(B0, V0, d, p_single);
    PenningTrap t_double = PenningTrap(B0, V0, d, p_double);


    //Simulations
    bool interact;
    interact = false;

    std::string outfile_s = "./res/out/single.bin";
    simulate(t_single, time, steps, outfile_s, interact, false);


    // interact = true;
    // std::string outfile_d = "./res/out/double.bin";
    // simulate(t_double, time, steps, outfile_d, interact, false);

    // Fill random 100
    // Vector of particles, double simulation
    // std::vector<Particle> p_random = fill_rand(100, d, q, m);
    // PenningTrap t_random = PenningTrap(B0, V0, d, p_random);



    // Relative error simulations
    std::string outfile;

    arma::mat a1 = analytical(r1,v1,B0,V0,d,m,q,4000.0,time);
    outfile = "./res/out/re_AN/4000.bin";
    a1.save(outfile, arma::arma_binary);
    outfile = "./res/out/re_RK4/4000.bin";
    simulate(t_single, time, 4000, outfile, false, false);
    outfile = "./res/out/re_FE/4000.bin";
    simulate(t_single, time, 4000, outfile, false, true);

    arma::mat a2 = analytical(r1,v1,B0,V0,d,m,q,8000.0,time);
    outfile = "./res/out/re_AN/8000.bin";
    a2.save(outfile, arma::arma_binary);
    outfile = "./res/out/re_RK4/8000.bin";
    simulate(t_single, time, 8000, outfile, false, false);
    outfile = "./res/out/re_FE/8000.bin";
    simulate(t_single, time, 8000, outfile, false, true);

    arma::mat a3 = analytical(r1,v1,B0,V0,d,m,q,16000.0,time);
    outfile = "./res/out/re_AN/16000.bin";
    a3.save(outfile, arma::arma_binary);
    outfile = "./res/out/re_RK4/16000.bin";
    simulate(t_single, time, 16000, outfile, false, false);
    outfile = "./res/out/re_FE/16000.bin";
    simulate(t_single, time, 16000, outfile, false, true);

    arma::mat a4 = analytical(r1,v1,B0,V0,d,m,q,32000.0,time);
    outfile = "./res/out/re_AN/32000.bin";
    a4.save(outfile, arma::arma_binary);
    outfile = "./res/out/re_RK4/32000.bin";
    simulate(t_single, time, 32000, outfile, false, false);
    outfile = "./res/out/re_FE/32000.bin";
    simulate(t_single, time, 32000, outfile, false, true);


    

 
    return 0;
}


void simulate(
    PenningTrap trap,
    double time,
    double steps,
    std::string outfile,
    bool interact,
    bool FE)
{
    // Timestep
    double dt = time/steps;

    // Out vectors

    // Out matrix


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




    // Initial values
    t_out(0) = 0.0;

    rx_out(0) = trap.p[0].r(0);
    ry_out(0) = trap.p[0].r(1);
    rz_out(0) = trap.p[0].r(2);

    vx_out(0) = trap.p[0].v(0);
    vy_out(0) = trap.p[0].v(1);
    vz_out(0) = trap.p[0].v(2);



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


    // Out filing



    // Write to file

    M_out.save(outfile, arma::arma_binary);


}

// Fill particle vector with random particles
std::vector<Particle> fill_rand(int n_particles, double d, double q, double m)
{
    std::vector<Particle> particles;

    for (int i = 0; i< n_particles; i++)
    {
        arma::vec r = arma::vec(3).randn() * 0.1 * d;  // random initial position
        arma::vec v = arma::vec(3).randn() * 0.1 * d;  // random initial velocity
        Particle p_i = Particle(q,m,r,v);

        particles.push_back(p_i);
        
    }
    return particles;
}