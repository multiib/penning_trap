#include <armadillo>
#include <stdexcept>
#include <cmath>

// Analtytical solution 
arma::mat analytical(
	arma::vec r_in, // Initial position (y value has to be 0)
	arma::vec v_in, // Initial velocity (x and z value has to be 0)
	double B0,      // Magnetic field strength
	double V0,      // Potential applied to the electrodes
	double d,       // Characteristic dimension
	double m,       // Mass of particle
	double q,       // Charge of particle
	double n,       // Number of steps
	double t)       // Time of simulation
{
	// Tests to check if initial values are legal
	if (r_in(1) != 0.0)
	{
		throw std::invalid_argument("Initial y position is not 0");
	}

	if (v_in(0) != 0.0)
	{
		throw std::invalid_argument("Initial x velocity is not 0.");
	}

	if (v_in(2) != 0.0)
	{
		throw std::invalid_argument("Initial z velocity is not 0.");
	}

	// Setup variables
	double dt = t/n;
	double x_0 = r_in(0);
	double z_0 = r_in(2);
	double v_0 = v_in(1);
	double w_0 = (q*B0)/m;
	double w_z = std::sqrt((2*q*V0)/(m*d*d));

	double w_plus  = (w_0 + std::sqrt(w_0*w_0 - 2*w_z*w_z))/2;
	double w_minus = (w_0 - std::sqrt(w_0*w_0 - 2*w_z*w_z))/2;

	// Amplitudes
	double A_plus = (v_0+w_minus*x_0)/(w_minus-w_plus);
	double A_minus = -(v_0+w_plus*x_0)/(w_minus-w_plus);

	// Out matrix
	arma::mat r = arma::mat(4, n);

	for (int i = 0; i < n; i++)
	{
		// Time vector
		r(0, i) = dt*i;

		// x vector
		r(1, i) = A_plus*std::cos(w_plus*dt*i)+A_minus*std::cos(w_minus*dt*i);

		// y vector
		r(2, i) = -(A_plus*std::sin(w_plus*dt*i)+A_minus*std::sin(w_minus*dt*i));

		//z vector
		r(3, i) = z_0*std::cos(w_z * dt*i);
	}

	return r;
}