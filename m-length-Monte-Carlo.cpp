#include "spa_library.h"   //contains the functions for main program.
#include "extra.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

ofstream outfile_mlength;
ofstream outfile_mcdetails;
ofstream placeholder;
ofstream outfile_freeenergy;

int main(int argc, char* argv[])
{
  int no_sweeps;
  double final_temp, initial_temp, temp_step;

  cout << "Enter the size, U: ";
  cin >> size >> U;
  cout << "Enter the number of MC sweeps, initial and final, temp_step: ";
  cin >> no_sweeps >> initial_temp >> final_temp >> temp_step;

  milliseconds begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  MatrixXcd Mc=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd Mcx=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd Mcy=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd Mcz=MatrixXcd::Zero(2*size,2*size);
  MatrixXd randsigma=MatrixXd::Zero(size,3);

  construct_h0(Mc);
  long idum = time(NULL);
  for(int i=0; i<size; i++) randsigma(i,2) = pow(-1,i); //sigma_generate(randsigma, i, idum);
  matrixelement_sigmax(Mcx, randsigma, size);
  matrixelement_sigmay(Mcy, randsigma, size);
  matrixelement_sigmaz(Mcz, randsigma, size);
  MatrixXcd initial_Hamiltonian = Mc-U/2*(Mcx+Mcy+Mcz);
  double free_energy = find_free_energy(initial_Hamiltonian, final_temp, randsigma, U, placeholder);

  MatrixXd suggested_randsigma = randsigma;
  MatrixXcd suggested_Mcy=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd suggested_Mcx=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd suggested_Mcz=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd suggested_Hamiltonian=MatrixXcd::Zero(2*size,2*size);
  double suggested_free_energy;

  string filename, latticedata;
  latticedata = "_U="+to_string(int(U))+"_size="+to_string(size);
  filename="data/m_length_vs_temp"+ current_time_str()+latticedata+".txt";  outfile_mlength.open(filename);
  filename="data/mcdetails"+current_time_str()+latticedata+".txt"; outfile_mcdetails.open(filename);
  filename="data/free_energy_vs_temp"+current_time_str()+latticedata+".txt"; outfile_freeenergy.open(filename);
  cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";

  assert(no_sweeps%2==0);
  double final_free_energy, count_free_energy;
  double uniform_rv, move_prob, m_length;

  for(double temperature=final_temp; temperature>=initial_temp; temperature-=temp_step)
  {
    final_free_energy = count_free_energy = 0;
    for(int sweep=0; sweep<no_sweeps; sweep++)
    {
      for(int lattice_index=0; lattice_index<size; lattice_index++)
      {
        sigma_generate(suggested_randsigma, lattice_index, idum);
        matrixelement_sigmax(suggested_Mcx, suggested_randsigma, size);
        matrixelement_sigmay(suggested_Mcy, suggested_randsigma, size);
        matrixelement_sigmaz(suggested_Mcz, suggested_randsigma, size);
        suggested_Hamiltonian = Mc-U/2*(suggested_Mcx+suggested_Mcy+suggested_Mcz);
        suggested_free_energy = find_free_energy(suggested_Hamiltonian, temperature, suggested_randsigma, U, placeholder);

        uniform_rv = ran0(&idum); move_prob = exp((free_energy - suggested_free_energy)/temperature);
        if( uniform_rv <= move_prob)
        {
          outfile_mcdetails << "accepted: Temp=" << temperature << " Sweep= " << sweep << " uniform_rv="<< uniform_rv << " move_prob=" << move_prob << " free_energy=" << free_energy << " suggested_free_energy= " << suggested_free_energy << endl;
          free_energy = suggested_free_energy;
          randsigma = suggested_randsigma;
        }
        else
        {
          outfile_mcdetails << "rejected: Temp=" << temperature << " Sweep= " << sweep << " uniform_rv="<< uniform_rv << " move_prob=" << move_prob << " free_energy=" << free_energy << " suggested_free_energy= " << suggested_free_energy << endl;
          suggested_randsigma=randsigma;
        }
      }
    }
    outfile_mlength << temperature << " ";
    double m_length_avg=0;
    for(int j=0; j<size; j++)
    {
      m_length = 0;
      for( int k=0; k<3; k++) m_length += pow(randsigma(j,k),2);
      outfile_mlength << " " << sqrt(m_length)  << " ";
      m_length_avg+= sqrt(m_length);
    }
    outfile_mlength << " \t" << m_length_avg/double(size) <<  endl;
    outfile_freeenergy << temperature << " " << free_energy << endl;
    progress_percent_desc(initial_temp, final_temp, temperature);
  }

  cout << endl;
  milliseconds end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  show_time(begin_ms, end_ms,"MC calculation");
  spinarrangement_Mathematica_output(randsigma);

  outfile_mcdetails.close();
  outfile_mlength.close();
  outfile_freeenergy.close();
  return 0;
}
