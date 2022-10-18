#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <armadillo>

#include "../lib/particle.hpp"
#include "../lib/penningtrap.hpp"

int main()
{

    double q =1;
    double m = 100.0;
    arma::vec r = arma::vec("20.0 0.0 20.0");
    arma::vec v = arma::vec("20.0 25.0 0.0");
    Particle first = Particle(q,m,r,v);


    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500;
    std::vector<Particle> p;
    p.push_back(first);

    PenningTrap trap = PenningTrap(B0, V0, d, p);

    // Setup simulation variables
    double time = 50.0; // microseconds
    int steps = 100;
    double dt = time/steps;


    // Out vectors
    arma::vec t_out = arma::vec(steps).fill(0.);
    arma::vec z_out = arma::vec(steps);


    
    z_out(0) = trap.p[0].r(2); //Initial z value

    for (int i = 1; i < steps; i++)
    {
        trap.evolve_forward_Euler(dt);
        // t_out(i) = i*dt;
        // z_out(i) = trap.p[0].r(2);
        // std::cout<<trap.p[0].r<<std::endl;
    }

    // arma::mat M_out = arma::mat(steps, 2);

    // M_out.col(0) = t_out;
    // M_out.col(1) = z_out;


    // //Write to file
    // std::string outfile = "./out/single.bin";
    // std::ofstream ofile;
    // M_out.save(outfile, arma::raw_ascii);




    return 0;
}