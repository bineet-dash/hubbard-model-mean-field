#include "tca_library.hpp"
#include "extra.hpp"

int main(int argc, char* argv[])
{

  cout << "Enter the size, U: ";
  cin >> size >> U;

  // int final_exp, initial_exp;
  // cout << "Enter the number of MC sweeps, final and initial exponent: ";
  // cin >> no_sweeps >> final_exp >> initial_exp;
  int no_sweeps = 20;
  int initial_exp = -2;
  int final_exp = -1;

  double final_temp = 10*pow(10,final_exp);

  milliseconds begin_ms, end_ms;
  begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  MatrixXd randsigma = MatrixXd::Zero(size, 3);
  long idum = time(NULL);
  for(int i=0; i<size; i++)  randsigma(i,2) = 1;
  for(int i=0; i<size; i++)  greens_sigma_generate(randsigma, i, idum);

  int Lc = size/4;

  string filename, latticedata;
  latticedata = "_U="+to_string(int(U))+"_size="+to_string(size)+"_sweeps="+to_string(no_sweeps)+"_Lc="+to_string(Lc);
  filename="wolframscripts/tca_spin_arrangement"+current_time_str()+latticedata+".nb"; ofstream outfile_spinarr(filename);
  spinarrangement_Mathematica_output(randsigma,outfile_spinarr);
  filename="data/tca_m_length_vs_temp"+ current_time_str()+latticedata+".txt"; ofstream outfile_mlength(filename);
  // filename="data/free_energy_vs_temp"+current_time_str()+latticedata+".txt"; ofstream outfile_freeenergy(filename);
  // filename="data/mcdetails"+current_time_str()+latticedata+".txt"; ofstream outfile_mcdetails(filename);
  cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";

  for(int j=final_exp; j>=initial_exp; j--)
  {
    for(double i=10; i>=2; i-=1)
    {
      double temperature = i*pow(10,j);
      for(int sweep=0; sweep<0.75*no_sweeps; sweep++)
      {
        for(int lattice_index=0; lattice_index<size; lattice_index++)
        {
          MatrixXd selected_randsigma = MatrixXd::Zero(Lc,3);
          for(int i=0; i<Lc; i++) selected_randsigma.row(i) = randsigma.row((lattice_index+i)%size);
          MatrixXcd H_cl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_randsigma);
          double free_energy = //find_canonical_free_energy(H_cl, temperature, randsigma);
                                filled_E(H_cl,Lc);

          MatrixXd suggested_randsigma_selected = selected_randsigma;
          greens_sigma_generate(suggested_randsigma_selected,0, idum); //present lattice index is 0 in selected part
          MatrixXcd suggested_H_cl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,suggested_randsigma_selected);
          double suggested_free_energy = //find_canonical_free_energy(suggested_H_cl, temperature, randsigma);
                                        filled_E(suggested_H_cl,Lc);

          double uniform_rv = ran0(&idum); double move_prob = exp((free_energy - suggested_free_energy)/temperature);
          if(uniform_rv <= move_prob)
          {
            free_energy = suggested_free_energy;
            randsigma.row(lattice_index) = suggested_randsigma_selected.row(0);
          }
          else
          {
            suggested_randsigma_selected=selected_randsigma;
          }
        }
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      for(int sweep= int(0.75*no_sweeps); sweep<no_sweeps; sweep++)
      {
        for(int lattice_index=0; lattice_index<size; lattice_index++)
        {
          MatrixXd selected_randsigma = MatrixXd::Zero(Lc,3);
          for(int i=0; i<Lc; i++) selected_randsigma.row(i) = randsigma.row((lattice_index+i)%size);
          MatrixXcd H_cl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_randsigma);
          double free_energy = //find_canonical_free_energy(H_cl, temperature, randsigma);
                               filled_E(H_cl,Lc);

          MatrixXd suggested_randsigma_selected = selected_randsigma;
          greens_sigma_generate(suggested_randsigma_selected,0, idum); //present lattice index is 0 in selected part
          MatrixXcd suggested_H_cl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,suggested_randsigma_selected);
          double suggested_free_energy =// find_canonical_free_energy(suggested_H_cl, temperature, randsigma);
                                        filled_E(suggested_H_cl,Lc);

          double uniform_rv = ran0(&idum); double move_prob = exp((free_energy - suggested_free_energy)/temperature);
          if(uniform_rv <= move_prob)
          {
            free_energy = suggested_free_energy;
            randsigma.row(lattice_index) = suggested_randsigma_selected.row(0);
          }
          else
          {
            suggested_randsigma_selected=selected_randsigma;
          }
        }
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      outfile_mlength << setw(5) << temperature << " ";
      for(int j=0; j<size; j++) outfile_mlength << " " << setw(5) << randsigma(j,2)  << " ";
      outfile_mlength << endl;
      // outfile_mlength << " \t" << m_length_avg/double(size) <<  endl;
      // outfile_freeenergy << temperature << " " << final_free_energy/double(count_free_energy) << endl;

      cout << "\rtemperature = " << temperature << " done."; cout.flush();
    }
  }

  cout << endl;
  end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  show_time(begin_ms, end_ms,"TCA calculation");
  spinarrangement_Mathematica_output(randsigma,outfile_spinarr);

  outfile_mlength.close();
  // outfile_mcdetails.close();
  // outfile_freeenergy.close();
  outfile_spinarr.close();
  return 0;
}
