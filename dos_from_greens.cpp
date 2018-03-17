#include "tca_library.hpp"
#include "extra.hpp"
#include <iomanip>

MatrixXcd cluster_h0_wo_pbc(int Lc)
{
  MatrixXcd Mc = MatrixXcd::Zero(2*Lc,2*Lc);
  for(int row=0; row<2*Lc-1; row++) Mc(row,row+1)=Mc(row+1,row)=-t;
  Mc(Lc,Lc-1)= Mc(Lc-1,Lc)=0;
  return Mc;
}

int main()
{
  std::cout << "Enter size and U: ";
  cin >> size >> U;

  MatrixXd randsigma = MatrixXd::Zero(size, 3);
  long idum = -1; //time(NULL);
  for(int i=0; i<size; i++)  randsigma(i,2) = 1;
  // for(int i=0; i<size; i++)  greens_sigma_generate(randsigma, i, idum);
  cout << "spins\n\n" << randsigma << endl << endl;

  MatrixXcd H = cluster_h0_wo_pbc(size)- U/2*cluster_sigmaz(size,randsigma);

  // MatrixXcd H_part = MatrixXcd::Zero(size,size);
  // H_part.block(0,0,size/2,size/2)=H.block(0,0,size/2,size/2);
  // H_part.block(size/2,size/2,size/2,size/2)=H.block(size,size,size/2,size/2);
  // cout << Eigenvalues(H_part).transpose() << endl << endl;
  // VectorXd eigenvals = Eigenvalues(H_part);


  ofstream fout("explicit_dos.txt");
  // ofstream fout2("explicit_dos_check.txt");
  for(double omega=omega_L; omega<omega_U; omega+=0.04)
  {
    MatrixXcd Z = cd(omega,eta)*MatrixXcd::Identity(H.rows(),H.rows());
    MatrixXcd G_total= invert(Z-H);
    MatrixXcd G_part = MatrixXcd::Zero(size,size);
    G_part.block(0,0,size/2,size/2)=G_total.block(0,0,size/2,size/2);
    G_part.block(size/2,size/2,size/2,size/2)=G_total.block(size,size,size/2,size/2);

    // MatrixXcd Z = cd(omega,eta)*MatrixXcd::Identity(H_part.rows(),H_part.rows());
    // MatrixXcd G_part = invert(Z-H_part);

    // cd sum = 0;
    // for(int i=0; i<eigenvals.size(); i++) sum+= cd(1,0)/(cd(omega,eta)-cd(eigenvals(i),0));
    // fout2 << omega << " " << -1/M_PI*sum.imag() << endl;

    double a_k_w = -1/M_PI*G_part.imag().trace();
    double a_k_w2 = -1/M_PI*G_total.imag().trace();
    fout << omega << " " << a_k_w << " " << a_k_w2 << endl;

    // cout.precision(2);
    // cout << G_total.imag() << endl << endl;
    // cout << G_part.imag() << endl << endl;
    // exit(1);
  }


}
