#ifndef GREEN_LIBRARY_HPP_INCLUDED
#define GREEN_LIBRARY_HPP_INCLUDED

#include "spa_library.hpp"
#include "extra.hpp"
#include <iomanip>

milliseconds begin_ms, end_ms;

double omega_L = -8;
double omega_U = 8;

double eta = 0.05;
double a = 1;


MatrixXcd invert(MatrixXcd A)
{
  int N = A.rows();
  int *IPIV = new int[N+1];
  int LWORK = N*N;
  __complex__ double* WORK= new __complex__ double [LWORK];
  int INFO;

  zgetrf_(&N,&N, reinterpret_cast <__complex__ double*> (A.data()), &N,IPIV,&INFO);
  zgetri_(&N, reinterpret_cast <__complex__ double*> (A.data()), &N,IPIV,WORK,&LWORK,&INFO);

  delete IPIV;
  delete WORK;
  if(INFO==0) return A;
  else {cout << "Inversion failed. Exiting..."; exit(28);}
}

double spectral_weight(MatrixXcd H, double k, double omega)
{
  MatrixXcd Z = cd(omega,eta)*MatrixXcd::Identity(H.rows(),H.rows());
  MatrixXcd G = invert(Z-H);
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

double check_e(MatrixXcd H)
{
  double ne=0.0;
  for(double omega=omega_L; omega<omega_U; omega += 0.01)
  {
    ne += dos(H,omega)*0.01;
  }
  cout << ne << endl;
  exit(1);
}

double filled_E(MatrixXcd H)
{
  double step = 0.02;
  double energy = 0.0;
  double no_of_electrons = 0.0;
  double omega;
  // begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  for(omega = omega_L; omega <= omega_U; omega += step)
  {
    if(no_of_electrons < size)
    {
      energy += omega*dos(H, omega)*step;
      no_of_electrons += dos(H, omega)*step;

      // end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
      // cout << setw(5) << omega << " " << setw(5) << (end_ms.count()-begin_ms.count()) << endl;
    }
    else break;
  }
  return energy;
}


#endif
