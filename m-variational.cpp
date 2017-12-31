#include "spa_library.h"
// #include "extra.h"

void find_dos(MatrixXcd Hamiltonian)
{
  ofstream fout("DoS.nb");
  int bin;
  cout << "Enter the bin-size:" << '\n';
  cin >> bin;

  ComplexEigenSolver <MatrixXcd> ces;
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


void check_tb_validity(MatrixXcd Hamiltonian, double t)
{
  ComplexEigenSolver <MatrixXcd> ces;
  ces.compute(Hamiltonian);

  cout << "The eigenvalues obtained numerically: \n";
  cout << ces.eigenvalues().real().transpose() << endl;

  int size=Hamiltonian.rows()/2;
  cout << "The eigenvalues given by theoretical TB model: \n";
  for(int i=0; i<size; i++) cout << -2*t*cos(2*M_PI*i/double(size)) << ", ";
  cout << endl;
}

double zero_temp_free_energy(MatrixXcd& H, MatrixXd sigma, ofstream& placeholder)
{
  double zero_temp_free_energy = U/4*sigma.unaryExpr(&Sqr).sum();
  std::vector<double> eigenvalues; diagonalize(H, eigenvalues);
  sort(eigenvalues.begin(),eigenvalues.end());
  for(int i=0; i<size/2; i++) zero_temp_free_energy += eigenvalues.at(i);
  return zero_temp_free_energy;
}

int main(int argc, char* argv[])
{
  cout << "Enter the size and U: ";
  cin >> size >> U;

  double min_m, max_m, step_m;
  cout << "Enter minimum m, maximum m, step size: ";
  cin >> min_m >> max_m >> step_m;

  // double Umin, Umax, Ustep;

  MatrixXcd Mc=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd Mcx=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd Mcy=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd Mcz=MatrixXcd::Zero(2*size,2*size);
  MatrixXd sigma=MatrixXd::Zero(size, 3);

  construct_h0(Mc);

  string filename, latticedata;
  latticedata = "_U="+to_string(int(U))+"_size="+to_string(size); ofstream placeholder;
  filename = "data/variational_data"+current_time_str()+latticedata+".txt"; ofstream outfile_variational(filename);

  // for(U=Umin; U<Umax; U+=Ustep)
  // {
    // for(int i=0; i<sigma.rows(); i++) sigma(i,2) = pow(-1,i)*min_m ;
    // matrixelement_sigmax(Mcx, sigma, size);
    // matrixelement_sigmay(Mcy, sigma, size);
    // matrixelement_sigmaz(Mcz, sigma, size);
    // MatrixXcd Hamiltonian = Mc-U/2*(Mcx+Mcy+Mcz);
    // free_energy = find_free_energy(Hamiltonian, temperature, sigma, U, placeholder);

    for(double m=min_m; m<max_m; m+= step_m)
    {
      for(int i=0; i<sigma.rows(); i++) sigma(i,2) = pow(-1,i)*m ;
      matrixelement_sigmax(Mcx, sigma, size);
      matrixelement_sigmay(Mcy, sigma, size);
      matrixelement_sigmaz(Mcz, sigma, size);
      MatrixXcd Hamiltonian = Mc-U/2*(Mcx+Mcy+Mcz);
      // double new_free_energy = zero_temp_free_energy(Hamiltonian, sigma, placeholder);
      double new_free_energy = find_internal_energy(Hamiltonian, sigma, placeholder);
      //if(new_free_energy<free_energy) {tmp=m; free_energy =new_free_energy;}
      outfile_variational << m << " " << new_free_energy << endl;
     }
    //progress_percent_asc(Umin,Umax,U);
  //}

  cout << endl;

  outfile_variational.close();
  return 0;
}
