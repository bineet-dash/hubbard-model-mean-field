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

string current_time_str(void)
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,sizeof(buffer),"%S-%M-%I-%Y-%m-%d",timeinfo);
  string str(buffer);
  return str;
}

double get_mu(double temperature, std::vector<double> v)
{
  sort (v.begin(), v.end());
  double bisection_up_lim = v.back();
  double bisection_low_lim = v.front();

  double mu, no_of_electrons; int count=0;
  double epsilon = 0.01;

  for(; ;)
  {
    int it=0; no_of_electrons=0;  count++;
    cout << "Loop:" << count << "\n----------------------\n"; if(count>10) exit(1);
    mu = 0.5*(bisection_low_lim+bisection_up_lim) ;

    while(no_of_electrons< double(size))
    {
      double fermi_func = 1/(exp((v[it]-mu)/temperature)+1); std::cout << fermi_func << '\t';
      no_of_electrons += fermi_func;  it++;
    }
    if(abs(no_of_electrons-size) < epsilon){cout << "exact " << v[it] << ", mu=" << mu << " position= "<< it << " " <<  no_of_electrons << endl; return mu; break;}
    else if(no_of_electrons > size+epsilon) { cout << "upper " << mu << " " << it-1 << " " <<  no_of_electrons << endl; bisection_up_lim=mu;}
    else if(no_of_electrons < size-epsilon) { cout << "lower " << mu << " " << it-1 << " " << no_of_electrons << endl; bisection_low_lim=mu;}
    cout << "\n-----------------------------\n";
  }
}

double get_mu(double temperature, VectorXd v)
{
  vector<double> stdv (v.data(),v.data()+v.size());
  return get_mu(temperature, stdv);
}

double debug_free_energy(MatrixXcd Mc, double temperature, MatrixXd randsigma, ofstream& debug)
{
  ComplexEigenSolver <MatrixXcd> ces;
  ces.compute(Mc);
  VectorXd sortedeivals=sortascending(ces.eigenvalues().real());

  cout << sortedeivals.transpose() << endl;

  double mu = get_mu(temperature,sortedeivals);
  cout << "temperature = " << temperature << ", mu = " << mu << endl;

  // debug << endl;
  // debug << "Eigenvalues :\n "<< setprecision(3) << sortedeivals.transpose() << endl;
  // debug << "progress: " << endl;

  double free_energy = 0; double ekt =0;

  for(int i=0; i<sortedeivals.size()/2; i++)
  {
    ekt = (sortedeivals(i)-mu)/temperature;
    if(!isinf(exp(-ekt))) free_energy += -temperature*log(1+exp(-ekt));
    else  free_energy += sortedeivals(i);
    // debug << free_energy << "-> ";
  }

  // debug << "free_energy (-kT lnZ) = " << free_energy << endl;
  free_energy += U/4*randsigma.unaryExpr(&Sqr).sum();
  // debug << "elastic cost= " << U/4*randsigma.unaryExpr(&Sqr).sum() << endl;
  // debug << " final_free_energy= " << free_energy << endl;
  // debug << endl;
  return free_energy;
}

double find_free_energy(MatrixXcd Mc, double temperature, MatrixXd randsigma)
{
  ComplexEigenSolver <MatrixXcd> ces;
  ces.compute(Mc);
  VectorXd sortedeivals=sortascending(ces.eigenvalues().real());

  double free_energy = 0; double ekt =0;

  for(int i=0; i<sortedeivals.size()/2; i++)
  {
    ekt = (sortedeivals(i))/temperature;
    if(!isinf(exp(-ekt))) free_energy += -temperature*log(1+exp(-ekt));
    else  free_energy += sortedeivals(i);
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

bool diagonalize(MatrixXcd A, vector<double>& lambda)
{
  int N = A.cols();
  if (A.rows()!=N)  return false;
  MatrixXcd v(N,N);

  // VectorXcd lambdac(N);
  int LDA = A.outerStride();
  int LDV = v.outerStride();
  int INFO = 0;
  cd* w = new cd [N];
  char Nchar = 'N';
  char Vchar = 'V';
  int LWORK = int(A.size())*4;
  VectorXcd WORK(LWORK);
  VectorXd RWORK(2*LDA);

  zgeev_(&Nchar, &Vchar, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, reinterpret_cast <__complex__ double*> (w), 0, &LDV, reinterpret_cast <__complex__ double*> (v.data()), &LDV,  reinterpret_cast <__complex__ double*> (WORK.data()), &LWORK, RWORK.data(), &INFO);

  // for(int i=0; i<N; i++) v.col(i)=v.col(i)/v.col(i).unaryExpr(&Sqr).sum();
  for(int i=0; i<N; i++) lambda.push_back(w[i].real());
  return INFO==0;
}


#endif
