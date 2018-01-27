#include "spa_library.hpp"   //contains the functions for main program.
#include "extra.hpp"
#include "green_library.hpp"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

void show_eigenvalues(MatrixXcd H)
{
  std::vector<double> eigenvalues;
  diagonalize(H, eigenvalues);
  for(auto it=eigenvalues.begin(); it!= eigenvalues.end(); it++) cout << *it << " ";
  cout << endl << endl;
}

void greens_sigma_generate(MatrixXd suggested_randsigma, int lattice_index, long & idum)
{
  if(ran0(&idum)<=0.5) suggested_randsigma(lattice_index,2) *= -1;
}


int main(int argc, char* argv[])
{

  cout << "Enter the size, U: ";
  cin >> size >> U;

  // int final_exp, initial_exp;
  // cout << "Enter the number of MC sweeps, final and initial exponent: ";
  // cin >> no_sweeps >> final_exp >> initial_exp;
  int no_sweeps = 500;
  int initial_exp = -2;
  int final_exp = 0;

  double final_temp = 10*pow(10,final_exp);

  milliseconds begin_ms, end_ms;

  MatrixXcd Mc,Mcx,Mcy,Mcz;
  Mc = Mcx = Mcy = Mcz = MatrixXcd::Zero(2*size,2*size);
  MatrixXd randsigma=MatrixXd::Zero(size,3);


  construct_h0(Mc);
  long idum = time(NULL);
  for(int i=0; i<size; i++)  randsigma(i,2) = 1;// pow(-1,i); //randsigma(i,2) = 5;

  for(int i=0; i<size; i++)  greens_sigma_generate(randsigma, i, idum);

  // matrixelement_sigmax(Mcx, randsigma);
  // matrixelement_sigmay(Mcy, randsigma);
  matrixelement_sigmaz(Mcz, randsigma);
  MatrixXcd initial_Hamiltonian = Mc-U/2*(Mcx+Mcy+Mcz);

  double free_energy = filled_E(initial_Hamiltonian);

  // show_eigenvalues(initial_Hamiltonian);
  // cout << "aux field:\n"  << randsigma << endl << endl;
  // double free_energy = find_free_energy(initial_Hamiltonian, 0.0001, randsigma);
  // cout << "energy = " << free_energy << endl;
  // exit(1);

  MatrixXd suggested_randsigma = randsigma;
  MatrixXcd suggested_Mc,suggested_Mcx,suggested_Mcy,suggested_Mcz;
  suggested_Mc = suggested_Mcx = suggested_Mcy = suggested_Mcz = MatrixXcd::Zero(2*size,2*size);

  string filename, latticedata;
  latticedata = "_U="+to_string(int(U))+"_size="+to_string(size)+"_sweeps="+to_string(no_sweeps);
  filename="wolframscripts/spin_arrangement"+current_time_str()+latticedata+".nb"; ofstream outfile_spinarr(filename);
  spinarrangement_Mathematica_output(randsigma,outfile_spinarr);
  filename="data/m_length_vs_temp"+ current_time_str()+latticedata+".txt"; ofstream outfile_mlength(filename);
  filename="data/free_energy_vs_temp"+current_time_str()+latticedata+".txt"; ofstream outfile_freeenergy(filename);
  // filename="data/mcdetails"+current_time_str()+latticedata+".txt"; ofstream outfile_mcdetails(filename);
  cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";

  for(int j=final_exp; j>=initial_exp; j--)
  {
    for(double i=10; i>=1; i-=1)
    {
      double temperature = i*pow(10,j);
      for(int sweep=0; sweep<0.75*no_sweeps; sweep++)
      {
        for(int lattice_index=0; lattice_index<size; lattice_index++)
        {
          greens_sigma_generate(suggested_randsigma,lattice_index, idum);
          matrixelement_sigmaz(suggested_Mcz, suggested_randsigma);
          MatrixXcd suggested_Hamiltonian = Mc-U/2*(suggested_Mcx+suggested_Mcy+suggested_Mcz);

          double suggested_free_energy = filled_E(suggested_Hamiltonian);

          double uniform_rv = ran0(&idum); double move_prob = exp((free_energy - suggested_free_energy)/temperature);
          if(uniform_rv <= move_prob)
          {
            free_energy = suggested_free_energy;
            randsigma = suggested_randsigma;
          }
          else
          {
            suggested_randsigma=randsigma;
          }
          exit(1);
        }
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      double final_free_energy = 0; double count_free_energy = 0;
      for(int sweep= int(0.75*no_sweeps); sweep<no_sweeps; sweep++)
      {
        for(int lattice_index=0; lattice_index<size; lattice_index++)
        {
          greens_sigma_generate(suggested_randsigma,lattice_index, idum);
          matrixelement_sigmax(suggested_Mcx, suggested_randsigma);
          matrixelement_sigmay(suggested_Mcy, suggested_randsigma);
          matrixelement_sigmaz(suggested_Mcz, suggested_randsigma);
          MatrixXcd suggested_Hamiltonian = Mc-U/2*(suggested_Mcx+suggested_Mcy+suggested_Mcz);
          double suggested_free_energy = filled_E(suggested_Hamiltonian);

          double uniform_rv = ran0(&idum); double move_prob = exp((free_energy - suggested_free_energy)/temperature);
          if(uniform_rv <= move_prob)
          {
            free_energy = suggested_free_energy;
            randsigma = suggested_randsigma;
          }
          else
          {
            suggested_randsigma=randsigma;
          }
        }
        if(sweep%50==0) {final_free_energy += free_energy; count_free_energy++;}
      }

      outfile_mlength << temperature << " ";
      double m_length_avg=0;
      for(int j=0; j<size; j++)
      {
        double m_length = 0;
        for( int k=0; k<3; k++) m_length += pow(randsigma(j,k),2);
        outfile_mlength << " " << sqrt(m_length)  << " ";
        m_length_avg+= sqrt(m_length);
      }
      outfile_mlength << " \t" << m_length_avg/double(size) <<  endl;
      outfile_freeenergy << temperature << " " << final_free_energy/double(count_free_energy) << endl;

      cout << "\rtemperature = " << temperature << " done."; cout.flush();
      // progress_percent_desc(initial_temp, final_temp, temperature);
    }
  }

  cout << endl;
  end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  show_time(begin_ms, end_ms,"MC calculation");
  spinarrangement_Mathematica_output(randsigma,outfile_spinarr);

  // outfile_mcdetails.close();
  outfile_mlength.close();
  outfile_freeenergy.close();
  outfile_spinarr.close();
  return 0;
}
