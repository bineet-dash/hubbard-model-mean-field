#include "spa_library.hpp"

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

double m_vs_f(double m, char c)
{
  MatrixXcd Mc,Mcx,Mcy,Mcz;
  Mc = Mcx = Mcy = Mcz = MatrixXcd::Zero(2*size,2*size);
  MatrixXd sigma=MatrixXd::Zero(size, 3);

  construct_h0(Mc);
  for(int i=0; i<sigma.rows(); i++) sigma(i,2) = pow(-1,i)*m ;
  matrixelement_sigmax(Mcx, sigma);
  matrixelement_sigmay(Mcy, sigma);
  matrixelement_sigmaz(Mcz, sigma);
  MatrixXcd Hamiltonian = Mc-U/2*(Mcx+Mcy+Mcz);

  ofstream placeholder;
  // double new_free_energy = zero_temp_free_energy(Hamiltonian, sigma, placeholder);
  double new_free_energy = (c !='f')? find_internal_energy(Hamiltonian, sigma, placeholder):find_free_energy(Hamiltonian,0.001,sigma,placeholder);
  return new_free_energy;
}

void plot(string filename)
{
  FILE *pipe = popen("gnuplot -p", "w");
  string s_cmd = "plot \""+ filename +"\" using 1:2 with lines \n";
  const char* cmd = s_cmd.c_str();
  fprintf(pipe,"%s", cmd);               // plot type
  fprintf(pipe, "exit \n");             // termination character
  fflush(pipe);                         // flush the pipe
  pclose(pipe);
}

int main(int argc, char* argv[])
{
  cout << "Enter the size and U: ";
  cin >> size >> U;

  double min_m, max_m, step_m;
  min_m = 0; max_m = 1.5; step_m = 0.01;
  // cout << "Enter minimum m, maximum m, step size: ";
  // cin >> min_m >> max_m >> step_m;

  string filename, latticedata;
  latticedata = "_U="+to_string(int(U))+"_size="+to_string(size);

  // for(; ;)
  // {
    filename = "data/variational_data"+current_time_str()+latticedata+".txt";
    ofstream outfile_variational(filename);
    // double temperature; cin >> temperature;

    double free_energy = m_vs_f(min_m,'i');
    for(double m=min_m; m<max_m; m+= step_m)
    {
      double new_free_energy = m_vs_f(m,'i');
      if(new_free_energy > free_energy) {cout << "minimum here: " << m << endl; break;}
      else {free_energy = new_free_energy;}
      outfile_variational << m << " " << new_free_energy << endl;
     }

    outfile_variational.close();
  // }




  return 0;
}
