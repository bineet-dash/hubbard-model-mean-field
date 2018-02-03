#ifndef _SPA_LIBRARY_HPP_INCLUDED_
#define _SPA_LIBRARY_HPP_INCLUDED_

#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <string>
#include <chrono>
#include <lapacke.h>
#include "common_globals.hpp"

using namespace std;
using namespace Eigen;
using namespace std::literals;
using namespace std::chrono;

typedef std::complex<double> cd;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
double ran0(long *idum)
{
   long  k;
   double ans;

   *idum ^= MASK;
   k = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   ans=AM*(*idum);
   *idum ^= MASK;
   return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

double t=1;
double U;
int size;

double Sqr(cd z){return norm(z);}
double filter(double x) {if(abs(x)<1e-7) return 0.0; else return x;}
void filter(std::vector<double>& v) {for(int i=0; i<v.size(); i++)  v[i]=filter(v[i]); }
void filter(VectorXd& v) {for(int i=0; i<v.size(); i++)  v(i)=filter(v(i));}

bool diagonalize(MatrixXcd A, vector<double>& lambda)
{
  int N = A.cols();
  int LDA = A.outerStride();
  int INFO = 0;
  double* w = new  double [N];
  char Nchar = 'N';
  char Vchar = 'V';
  char Uchar = 'U';
  int LWORK = int(A.size())*4;
  __complex__ double* WORK= new __complex__ double [LWORK];
  double* RWORK = new double [3*LDA];
  zheev_( &Nchar, &Uchar, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, w, WORK, &LWORK, RWORK, &INFO );

  lambda.clear();
  for(int i=0; i<N; i++) lambda.push_back(w[i]);

  delete[] w; delete[] RWORK; delete[] WORK;
  return INFO==0;
}

void construct_h0(MatrixXcd &Mc)
{
  for(int row=0;row<2*size-1; row++) Mc(row,row+1)=Mc(row+1,row)=-t;
  Mc(size-1,0)=Mc(0,size-1)=-t; //PBC
  Mc(2*size-1,size)=Mc(size,2*size-1)= -t; //PBC
  Mc(size,size-1)= Mc(size-1,size)=0;
}

void sigma_generate(MatrixXd& randsigma, int i, long & idum, double T)
{
  double radius, u, theta;
  radius = (T<0.5)? 0.5+ran0(&idum): ( min(2*T,8.0)*ran0(&idum) );
  u = 2*ran0(&idum)-1;
  theta  = 2*3.1416*ran0(&idum);
  randsigma(i,0)= radius*sqrt(1-pow(u,2))*cos(theta);
  randsigma(i,1)= radius*sqrt(1-pow(u,2))*sin(theta);
  randsigma(i,2)= radius*u;
}

void matrixelement_sigmax(MatrixXcd &Mcx, MatrixXd randsigma)
{
  for(int row=0; row<size; row++)
  {
    Mcx(row,row+size) = cd(randsigma(row,0),0);
    Mcx(row+size,row) = cd(randsigma(row,0),0);
  }
}

void matrixelement_sigmay(MatrixXcd &Mcy, MatrixXd randsigma)
{
  for(int row=0; row<size; row++)
  {
    Mcy(row,row+size) = cd(0,-randsigma(row,1));
    Mcy(row+size,row) = cd(0,randsigma(row,1));
  }
}

void matrixelement_sigmaz(MatrixXcd &Mcz, MatrixXd randsigma)
{
  for(int row=0; row<size; row++)
    Mcz(row,row)= cd(randsigma(row,2),0);

  for(int row=size; row<2*size; row++)
    Mcz(row, row)=cd(-randsigma(row-size,2),0);
}

VectorXd sortascending(VectorXd v1)
{
  vector<double> stdv1 (v1.data(),v1.data()+v1.size());
  std::sort (stdv1.begin(), stdv1.end());
  Map<ArrayXd> sorted(stdv1.data(), stdv1.size());
  return sorted;
}

double get_mu(double temperature, std::vector<double> v)
{
  sort (v.begin(), v.end());
  double bisection_up_lim = v.back();
  double bisection_low_lim = v.front();

  double mu, no_of_electrons; int count=0;
  double epsilon = 0.000001;

  for(; ;)
  {
    no_of_electrons=0;
    mu = 0.5*(bisection_low_lim+bisection_up_lim);

    for(auto it = v.begin(); it!= v.end(); it++)
    {
      double fermi_func = 1/(exp((*it-mu)/temperature)+1);
      no_of_electrons += fermi_func;
    }
    if(abs(no_of_electrons-size) < epsilon)
    {
      return mu; break;
    }
    else if(no_of_electrons > size+epsilon)
    {
       if(abs(bisection_up_lim-v.front())<0.001){return mu; break;}
       else {bisection_up_lim=mu;}
    }
    else if(no_of_electrons < size-epsilon)
    {bisection_low_lim=mu;}
  }
}


double get_mu(double temperature, VectorXd v)
{
  vector<double> stdv (v.data(),v.data()+v.size());
  return get_mu(temperature, stdv);
}

double find_free_energy(MatrixXcd Mc, double temperature, MatrixXd randsigma)
{
  std::vector<double> eigenvalues;
  diagonalize(Mc, eigenvalues);
  sort(eigenvalues.begin(),eigenvalues.end());

  double free_energy = 0; double ekt =0;
  double mu = get_mu(temperature, eigenvalues);

  for(auto it=eigenvalues.begin(); it!= eigenvalues.end(); it++)
  {
    ekt = (*it-mu)/temperature;
    if(!isinf(exp(-ekt))) free_energy += -temperature*log(1+exp(-ekt));
    else  free_energy += (*it-mu);
  }

  free_energy += U/4*randsigma.unaryExpr(&Sqr).sum();
  return free_energy;
}

double find_canonical_free_energy(MatrixXcd Mc, double temperature, MatrixXd randsigma)
{
  std::vector<double> eigenvalues;
  diagonalize(Mc, eigenvalues);
  sort(eigenvalues.begin(),eigenvalues.end());

  double free_energy = 0; double ekt =0;

  for(auto it=eigenvalues.begin(); it!= eigenvalues.end(); it++)
  {
    ekt = (*it)/temperature;
    if(!isinf(exp(-ekt))) free_energy += -temperature*log(1+exp(-ekt));
    else  free_energy += (*it);
  }

  free_energy += U/4*randsigma.unaryExpr(&Sqr).sum();
  return free_energy;
}

double find_internal_energy(MatrixXcd Mc, MatrixXd randsigma, ofstream& debug)
{
  ComplexEigenSolver <MatrixXcd> ces;
  ces.compute(Mc);
  VectorXd sortedeivals=sortascending(ces.eigenvalues().real());

  double internal_energy = 0;
  for(int i=0; i<sortedeivals.size()/2; i++) internal_energy += sortedeivals(i);
  internal_energy += U/4*randsigma.unaryExpr(&Sqr).sum();
  return internal_energy;
}



#endif
