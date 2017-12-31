#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include "spalib.h"
#include "debug.h"

using namespace std;
using namespace Eigen;

void find_dos(MatrixXcf Hamiltonian)
{
  ofstream fout("DoS.nb");
  int bin;
  cout << "Enter the bin-size: ";
  cin >> bin;

  ComplexEigenSolver <MatrixXcf> ces;
  ces.compute(Hamiltonian);

  fout << "Histogram[{ ";
  for(int i=0; i<Hamiltonian.rows(); i++)
  {
    fout << ces.eigenvalues()[i].real();
    if(i!=Hamiltonian.rows()-1) fout << ", ";
  }
  fout << "}, " << bin << "]";
  fout.close();
}


void check_tb_validity(MatrixXcf Hamiltonian, float t)
{
  ComplexEigenSolver <MatrixXcf> ces;
  ces.compute(Hamiltonian);

  cout << "The eigenvalues obtained numerically: \n";
  cout << ces.eigenvalues().real().transpose() << endl;

  int size=Hamiltonian.rows()/2;
  cout << "The eigenvalues obtained by theoretical TB model: \n";
  for(int i=0; i<size; i++) cout << -2*t*cos(2*M_PI*i/float(size)) << ", ";
  cout << endl;
}

int main(int argc, char* argv[])
{
  int size, bin;
  float t, U, free_energy, internal_energy;

  cout << "Enter the size, t, U: ";
  cin >> size >> t  >> U ;


  MatrixXcf Mc=MatrixXcf::Zero(2*size,2*size);
  MatrixXcf Mcx=MatrixXcf::Zero(2*size,2*size);
  MatrixXcf Mcy=MatrixXcf::Zero(2*size,2*size);
  MatrixXcf Mcz=MatrixXcf::Zero(2*size,2*size);
  MatrixXf randsigma=MatrixXf::Zero(size, 3);

  construct_h0(Mc, t);

  for(int i=0; i<randsigma.rows(); i++) randsigma(i,2) = 1 ;

  matrixelement_sigmax(Mcx, randsigma, size);
  matrixelement_sigmay(Mcy, randsigma, size);
  matrixelement_sigmaz(Mcz, randsigma, size);
  MatrixXcf Hamiltonian = Mc-U/2*(Mcx+Mcy+Mcz);

  find_dos(Hamiltonian);

  //cout << Mc.real() << endl;

}
