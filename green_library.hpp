#ifndef GREEN_LIBRARY_HPP_INCLUDED
#define GREEN_LIBRARY_HPP_INCLUDED

#include "spa_library.hpp"

double omega_L = -8;
double omega_U = 8;

double eta = 0.05;
double a = 1;

// void get_Z( MatrixXcd& Z, double omega)
// {
// }

void greens(MatrixXcd H, MatrixXcd& G, double omega)
{
  MatrixXcd Z = cd(omega,eta)*MatrixXcd::Identity(H.rows(),H.rows());
  G = (Z-H).inverse();
}

double spectral_weight(MatrixXcd H, double k, double omega)
{
    MatrixXcd G;  greens(H, G,omega);
    double weight = G.trace().imag();
    return -1/M_PI*weight;
}

double dos(MatrixXcd H, double omega)
{
  double DoS = 0;
  for(int n=0; n<size; n++)
  {
    double k = n*M_PI/(size*a);// -M_PI/a ;
    DoS += spectral_weight(H, k, omega);
  }
  return DoS/size;
}

double mu(MatrixXcd H)
{
  double step = 0.02;
  double trapez_sum = 0.0; double omega = omega_L;
  while(trapez_sum < size)
  {
    trapez_sum += dos(H, omega)*step;
    omega += step;
  }
  return omega;
}

double filled_E(MatrixXcd H)
{
  double ef = mu(H);
  double step = 0.02;
  double energy = 0.0; double omega = 0;
  for(omega = omega_L; omega <= ef; omega += step)
  {
    energy += omega*dos(H, omega)*step;
    omega += step;
  }
  return energy;
}


#endif
