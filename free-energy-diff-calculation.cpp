#include "tca_library.hpp"
#include "extra.hpp"

int main(int argc, char* argv[])
{

  cout << "Enter the size, U: ";
  cin >> size >> U;

  // int final_exp, initial_exp;
  // cout << "Enter the number of MC sweeps, final and initial exponent: ";
  // cin >> no_sweeps >> final_exp >> initial_exp;
  int no_sweeps = 200;
  int initial_exp = -2;
  int final_exp = -1;

  double final_temp = 10*pow(10,final_exp);

  milliseconds begin_ms, end_ms;
  begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  MatrixXd randsigma = MatrixXd::Zero(size, 3);
  long idum = -1; //time(NULL);
  for(int i=0; i<size; i++)  randsigma(i,2) = 1;
  for(int i=0; i<size; i++)  greens_sigma_generate(randsigma, i, idum);

  int Lc = size/2;

  string filename, latticedata;
  latticedata = "_U="+to_string(int(U))+"_size="+to_string(size)+"_sweeps="+to_string(no_sweeps)+"_Lc="+to_string(Lc);
  filename="data/free_energy_diff_vs_temp_"+current_time_str()+latticedata+".txt"; ofstream outfile_freeenergy(filename);
  cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";

  for(int j=final_exp; j>=initial_exp; j--)
  {
    for(double i=10; i>=2; i-=1)
    {
      int count_free_energy = 0;
      double final_free_energy = 0;
      double temperature = i*pow(10,j);

      // // for(int lattice_index=0; lattice_index<size; lattice_index++)
      // // {
        int lattice_index = 2;

        // cout << randsigma.col(2).transpose() << endl;
        MatrixXd suggested_randsigma = randsigma;
        suggested_randsigma(lattice_index,2) *= -1; //flip the field moment
        // cout << suggested_randsigma.col(2).transpose() << endl;

        VectorXd free_energy_arr = VectorXd::Zero(4);
        free_energy_arr << filled_E_ed_disconnected(randsigma,temperature,lattice_index,Lc),
                           filled_E_ed_connected(randsigma,temperature,lattice_index,Lc),
                           filled_E_disconnected(randsigma, lattice_index, Lc, temperature),
                           filled_E_connected(randsigma, lattice_index, Lc, temperature);

        VectorXd changed_free_energy_arr = VectorXd::Zero(4);
        changed_free_energy_arr << filled_E_ed_disconnected(suggested_randsigma,temperature,lattice_index,Lc),
                                  filled_E_ed_connected(suggested_randsigma,temperature,lattice_index,Lc),
                                  filled_E_disconnected(suggested_randsigma, lattice_index, Lc, temperature),
                                  filled_E_connected(suggested_randsigma, lattice_index, Lc, temperature);

        VectorXd diff_free_energy_arr = changed_free_energy_arr - free_energy_arr;
      // }      
     
      outfile_freeenergy << temperature << " " << diff_free_energy_arr.transpose() << endl;
      cout << "\rtemperature = " << temperature << " done."; cout.flush();
    }
  }

  cout << endl;
  end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  show_time(begin_ms, end_ms,"TCA calculation");

  outfile_freeenergy.close();
  return 0;
}