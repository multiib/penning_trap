#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <armadillo>

#include "../lib/particle.hpp"
#include "../lib/penningtrap.hpp"






int main()
{

    // Simulation configurations
    double time = 50.0; // µs
    int steps = 100;



    // Penning Trap properties
    double B0 = 9.65e1; // 1.00 T = u/((µs)**2*e)
    double V0 = 2.41e6; // 25.o mV = u(µm)**2/((µs)**2*e)
    double d = 500; // µm


    // Ca+ properties
    double q = 1.0;    // Charge +1 e
    double m = 40.078; // Mass of Ca+


    // Particle 1
    arma::vec r = arma::vec("20.0 0.0 20.0"); // µm
    arma::vec v = arma::vec("0.0 25.0 0.0"); // µm/µs
    Particle p1 = Particle(q,m,r,v);


    // Particle 2
    arma::vec r = arma::vec("25.0 25.0 0.0"); // µm
    arma::vec v = arma::vec("0.0 40.0 5.0"); // µm/µs
    Particle p2 = Particle(q,m,r,v);


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


    // Simulations
    std::string outfile = "./res/out/single.bin";
    simulate(t_single, time, steps, outfile);

    std::string outfile = "./res/out/double.bin";
    simulate(t_single, time, steps, outfile);

 
    return 0;
}


void simulate(PenningTrap trap, double time, double steps, std::string outfile)
{
    // Timestep
    double dt = time/steps;

    // Out vectors
    // Time vec
    arma::vec t_out = arma::vec(steps);

    // Position vec
    arma::vec rx_out = arma::vec(steps);
    arma::vec ry_out = arma::vec(steps); 
    arma::vec rz_out = arma::vec(steps);

    // Velocity vec
    arma::vec vx_out = arma::vec(steps); 
    arma::vec vy_out = arma::vec(steps); 
    arma::vec vz_out = arma::vec(steps); 


    // Initial values
    t_out(0) = 0.0;

    rx_out(0) = trap.p[0].r(0);
    ry_out(0) = trap.p[0].r(1);
    rz_out(0) = trap.p[0].r(2);

    vx_out(0) = trap.p[0].v(0);
    vy_out(0) = trap.p[0].v(1);
    vz_out(0) = trap.p[0].v(2);


    // Running simulation
    for (int i = 1; i < steps; i++)
    {
        t_out(i) = i*dt;

        trap.evolve_RK4(dt);
        rx_out(i) = trap.p[0].r(0);
        ry_out(i) = trap.p[0].r(1);
        rz_out(i) = trap.p[0].r(2);
        vx_out(i) = trap.p[0].v(0);
        vy_out(i) = trap.p[0].v(1);
        vz_out(i) = trap.p[0].v(2);
    }


    // Out matrix
    arma::mat M_out = arma::mat(steps, 7);

    M_out.col(0) = t_out;
    M_out.col(1) = rx_out;
    M_out.col(2) = ry_out;
    M_out.col(3) = rz_out;
    M_out.col(4) = vx_out;
    M_out.col(5) = vy_out;
    M_out.col(6) = vz_out;

    // Write to file
    M_out.save(outfile, arma::raw_ascii);
}

//     