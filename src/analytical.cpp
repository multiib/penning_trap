#include <armadillo>
#include <stdexcept>
#include <cmath>

arma::mat analytical(
	arma::vec r_in,
	arma::vec v_in,
	double B0,
	double V0,
	double d,
	double m,
	double q,
	double n,
	double t)
{
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

	double dt = t/n;

	double x_0 = r_in(0);
	double z_0 = r_in(2);

	double v_0 = v_in(1);

	double w_0 = (q*B0)/m;
	double w_z = std::sqrt((2*q*V0)/(m*d*d));



	double w_plus  = (w_0 + std::sqrt(w_0*w_0 - 2*w_z*w_z))/2;
	double w_minus = (w_0 - std::sqrt(w_0*w_0 - 2*w_z*w_z))/2;

	double A_plus = (v_0+w_minus*x_0)/(w_minus-w_plus);
	double A_minus = -(v_0+w_plus*x_0)/(w_minus-w_plus);

	arma::mat r = arma::mat(4, n);

	for (int i = 0; i < n; i++)
	{
		r(0, i) = dt*i;
		r(1, i) = A_plus*std::cos(w_plus*dt*i)+A_minus*std::cos(w_minus*dt*i);
		r(2, i) = -(A_plus*std::sin(w_plus*dt*i)+A_minus*std::sin(w_minus*dt*i));
		r(3, i) = z_0*std::cos(w_z * dt*i);
	}


	return r;
}