#ifndef TCA_LIBARARY_HPP_INCLUDED
#define TCA_LIBARARY_HPP_INCLUDED

#include "green_library.hpp"

MatrixXcd cluster_h0(int Lc, int pbc_pref = 1)
{
  MatrixXcd Mc = MatrixXcd::Zero(2*Lc,2*Lc);
  for(int row=0; row<2*Lc-1; row++) Mc(row,row+1)=Mc(row+1,row)=-t;
  Mc(Lc,Lc-1)= Mc(Lc-1,Lc)=0; //up and down spins don't mix

  if (pbc_pref==1) {
    Mc(2*Lc-1,Lc)=Mc(Lc,2*Lc-1)= -t; //PBC
    Mc(Lc-1,0)=Mc(0,Lc-1)=-t; //PBC
  }
  return Mc;
}

MatrixXcd cluster_sigmaz(int Lc, MatrixXd rs)
{
  MatrixXcd Mcz = MatrixXcd::Zero(2*Lc,2*Lc);
  for(int row=0; row<Lc; row++)  Mcz(row,row)= cd(rs(row,2),0);
  for(int row=Lc; row<2*Lc; row++) Mcz(row, row)=cd(-rs(row-Lc,2),0);
  return Mcz;
}

void greens_sigma_generate(MatrixXd& suggested_randsigma, int lattice_index, long & idum)
{
  if(ran0(&idum)<=0.5) suggested_randsigma(lattice_index,2) *= -1;
}


#endif
