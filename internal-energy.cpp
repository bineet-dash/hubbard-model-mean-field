#include "spa_library.h"   //contains the functions for main program.
#include "extra.h"    //contains debugging functions (for outputting matrices and generaing Mathematica scripts)

int main(int argc, char* argv[])
{
  int size, no_sweeps;
  double t, U, free_energy, final_temp, initial_temp,  internal_energy;

  cout << "Enter the size, U: ";
  cin >> size >> U;
  cout << "Enter the number of MC sweeps, initial and final temperature: ";
  cin >> no_sweeps >> initial_temp >> final_temp;

/* ========================================================================================================================
  H_SPA= \sum_i{-t(c_j\dagger c_i+c_i\dagger c_j)+U/2*[n_i*<n_i>-\vec{m}_i.<\vec{m}_i>]+ U/4[<\vec{m}_i>^2-<n_i>^2]}-\mu N
  For half-filling, \mu=U/2.
  For half-filling and translationally invariant system, <n_i>=1 for each i
  So, H_SPA= \sum_i{-t(c_j\dagger c_i+c_i\dagger c_j)+ 3*U*N/4 - U*\vec{m}_i.<\vec{m}_i>/2 - U*<\vec{m}_i>^2/4
  We choose the spin quantization axis along z-axis.
  So H_SPA=H_0 + H_sample + H_elastic
  where, H_0= \sum_i{-t(c_j\dagger c_i+c_i\dagger c_j)+ 3*U*N/4
        H_sample = -U/2*(m_x*<m_x>+m_y<m_y>+m_z<m_z>)
        H_elastic = - U/4*<\vec{m}_i>^2
  We are going to diagonalize H_0+H_sample.
 =========================================================================================================================== */


  string filename, latticedata;
  latticedata = "_U="+to_string(int(U))+"_size="+to_string(size); ofstream placeholder;
  filename="data/mc_internal_energy_details"+current_time_str()+latticedata+".txt"; ofstream outfile_mcdetails(filename);
  filename="data/internal_energy_vs_temp"+current_time_str()+latticedata+".txt"; ofstream outfile_internalenergy(filename);
  cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";

  MatrixXcd Mc=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd Mcx=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd Mcy=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd Mcz=MatrixXcd::Zero(2*size,2*size);
  MatrixXd randsigma(size, 3);

  construct_h0(Mc); //Matrix for H_0
  long idum = time(NULL);
  for(int i=0; i<size; i++) sigma_generate(randsigma,i, idum);
  matrixelement_sigmax(Mcx, randsigma, size);
  matrixelement_sigmay(Mcy, randsigma, size);
  matrixelement_sigmaz(Mcz, randsigma, size);
  MatrixXcd initial_M = Mc-U/2*(Mcx+Mcy+Mcz);
  internal_energy = find_internal_energy(initial_M, final_temp, randsigma, U, placeholder);

  MatrixXd suggested_randsigma(size, 3);
  MatrixXcd suggested_Mcy=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd suggested_Mcx=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd suggested_Mcz=MatrixXcd::Zero(2*size,2*size);
  MatrixXcd suggested_Mc=MatrixXcd::Zero(2*size,2*size);
  double suggested_internal_energy;

  assert(no_sweeps%2==0);
  double final_internal_energy, count_internal_energy;
  double uniform_rv, move_prob;

  for(double temperature=final_temp; temperature>=initial_temp; temperature-=0.01)
  {
    final_internal_energy = count_internal_energy = 0;
    for(int i=0; i<no_sweeps; i++)
    {
      for(int lattice_index=0; lattice_index<size; lattice_index++)
      {
        sigma_generate(suggested_randsigma,lattice_index, idum);
        matrixelement_sigmax(suggested_Mcx, suggested_randsigma, size);
        matrixelement_sigmay(suggested_Mcy, suggested_randsigma, size);
        matrixelement_sigmaz(suggested_Mcz, suggested_randsigma, size);
        suggested_Mc = Mc-U/2*(suggested_Mcx+suggested_Mcy+suggested_Mcz);
        suggested_internal_energy = find_internal_energy(suggested_Mc, final_temp, suggested_randsigma, U, placeholder);

        uniform_rv = ran0(&idum); move_prob = exp((internal_energy - suggested_internal_energy)/temperature);
        if( uniform_rv <= move_prob)
        {
          outfile_mcdetails << "accepted: ran0=" << uniform_rv << " ; prob=" << move_prob << "; internal_energy=" << internal_energy << "; suggested_internal_energy= " << suggested_internal_energy << endl;
          internal_energy = suggested_internal_energy;
          randsigma = suggested_randsigma;
        }
        else
        {
          outfile_mcdetails << "rejected: ran0=" << uniform_rv << " ; prob=" << move_prob << "; internal_energy=" << internal_energy << "; suggested_internal_energy= " << suggested_internal_energy << endl;
        }
      }
      progress_percent_desc(initial_temp, final_temp, temperature);
    }
    cout << endl;
  }

  return 0;

}
