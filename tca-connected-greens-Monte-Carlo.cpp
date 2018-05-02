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
  long idum = time(NULL);
  for(int i=0; i<size; i++)  randsigma(i,2) = 1;
  for(int i=0; i<size; i++)  greens_sigma_generate(randsigma, i, idum);

  int Lc = size/2;

  string filename, latticedata;
  latticedata = "_U="+to_string(int(U))+"_size="+to_string(size)+"_sweeps="+to_string(no_sweeps)+"_Lc="+to_string(Lc);
  filename="wolframscripts/tca_connected_spin_arrangement"+current_time_str()+latticedata+".nb"; ofstream outfile_spinarr(filename);
  spinarrangement_Mathematica_output(randsigma,outfile_spinarr);
  filename="data/tca_connected_m_length_vs_temp"+ current_time_str()+latticedata+".txt"; ofstream outfile_mlength(filename);
  filename="data/tca_connected_free_energy_vs_temp"+current_time_str()+latticedata+".txt"; ofstream outfile_freeenergy(filename);
  // filename="data/tca_connected_dos"+current_time_str()+latticedata+".txt"; ofstream outfile_dos(filename);
  filename="data/mcdetails"+current_time_str()+latticedata+".txt"; ofstream outfile_mcdetails(filename);
  cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";

  for(int j=final_exp; j>=initial_exp; j--)
  {
    for(double i=10; i>=2; i-=1)
    {
      int count_free_energy = 0;
      double final_free_energy = 0;
      double temperature = i*pow(10,j);
      // double temperature = 0.001;
    
      // MatrixXd dos_l, dos_r;
      // MatrixXd average_dos_l = MatrixXd::Zero(dos_l.rows(),2);
      // MatrixXd average_dos_r = MatrixXd::Zero(dos_r.rows(),2);
      // MatrixXd updated_dos_l= MatrixXd::Zero(dos_l.rows(),2);
      // MatrixXd updated_dos_r= MatrixXd::Zero(dos_r.rows(),2);
      // int count_for_dos = 0;

      for(int sweep=0; sweep<0.75*no_sweeps; sweep++)
      {
        for(int lattice_index=0; lattice_index<size; lattice_index++)
        {
          // double free_energy = filled_E_connected(randsigma, lattice_index, Lc, dos_l, dos_r);
          double free_energy = filled_E_ed(randsigma,temperature,lattice_index,Lc);

          MatrixXd suggested_randsigma = randsigma;
          greens_sigma_generate(suggested_randsigma,lattice_index, idum); 
          //  cout << suggested_randsigma << endl;
          // double suggested_free_energy = filled_E_connected(suggested_randsigma, lattice_index, Lc, dos_l, dos_r);
          double suggested_free_energy = filled_E_ed(suggested_randsigma,temperature,lattice_index,Lc);

          double uniform_rv = ran0(&idum); double move_prob = exp((free_energy - suggested_free_energy)/temperature);
          if(uniform_rv <= move_prob)
          {
            free_energy = suggested_free_energy;
            randsigma.row(lattice_index) = suggested_randsigma.row(lattice_index);
            outfile_mcdetails << "accepted: Temp=" << temperature << " Sweep= " << sweep << " uniform_rv="<< uniform_rv << " move_prob=" << move_prob << " free_energy=" << free_energy << " suggested_free_energy= " << suggested_free_energy << endl; 
          }
          else
          {
            outfile_mcdetails << "rejected: Temp=" << temperature << " Sweep= " << sweep << " uniform_rv="<< uniform_rv << " move_prob=" << move_prob << " free_energy=" << free_energy << " suggested_free_energy= " << suggested_free_energy << endl;
          }
        }
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      for(int sweep= int(0.75*no_sweeps); sweep<no_sweeps; sweep++)
      {
        for(int lattice_index=0; lattice_index<size; lattice_index++)
        {
          // double free_energy = filled_E_connected(randsigma, lattice_index, Lc, dos_l, dos_r);
          double free_energy = filled_E_ed(randsigma,temperature,lattice_index,Lc);

          MatrixXd suggested_randsigma = randsigma;
          greens_sigma_generate(suggested_randsigma,lattice_index, idum); 
          // cout << suggested_randsigma << endl;
          // double suggested_free_energy = filled_E_connected(suggested_randsigma, lattice_index, Lc, updated_dos_l, updated_dos_r);
          double suggested_free_energy = filled_E_ed(suggested_randsigma,temperature,lattice_index,Lc);

          double uniform_rv = ran0(&idum); double move_prob = exp((free_energy - suggested_free_energy)/temperature);

          if(uniform_rv <= move_prob)
          {
            free_energy = suggested_free_energy;
            randsigma.row(lattice_index) = suggested_randsigma.row(lattice_index);
            outfile_mcdetails << "accepted: Temp=" << temperature << " Sweep= " << sweep << " uniform_rv="<< uniform_rv << " move_prob=" << move_prob << " free_energy=" << free_energy << " suggested_free_energy= " << suggested_free_energy << endl;
            
            if(sweep%5==0)
            {
              final_free_energy += free_energy;
            }
          }
          else
          {
            final_free_energy += free_energy;
            outfile_mcdetails << "rejected: Temp=" << temperature << " Sweep= " << sweep << " uniform_rv="<< uniform_rv << " move_prob=" << move_prob << " free_energy=" << free_energy << " suggested_free_energy= " << suggested_free_energy << endl;
          }

          // if(uniform_rv <= move_prob){
          //   free_energy = suggested_free_energy;
          //   randsigma.row(lattice_index) = suggested_randsigma.row(lattice_index);
          //   if(sweep%5==0){
          //     average_dos_l += updated_dos_l;
          //     average_dos_r += updated_dos_r;
          //     count_for_dos++;
          //   }
          // }
          // else{
          //   if(sweep%5==0){ 
          //     average_dos_l += updated_dos_l;
          //     average_dos_r += updated_dos_r;
          //     count_for_dos++;
          //   }
          // }
        }
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      outfile_mlength << setw(5) << temperature << " ";
      for(int j=0; j<size; j++) outfile_mlength << " " << setw(5) << randsigma(j,2)  << " ";
      outfile_mlength << endl;

      outfile_freeenergy << temperature << " " << final_free_energy/double(count_free_energy) << endl;
      // for(int i=0; i<average_dos_l.rows(); i++)
      // {
      //   outfile_dos << average_dos_l(i,0)/double(count_for_dos) << " " << average_dos_l(i,1)/double(count_for_dos) << " " <<  average_dos_r(i,1)/double(count_for_dos) << endl;
      // }
      // cout << "count for dos = " << count_for_dos << endl;

      cout << "\rtemperature = " << temperature << " done."; cout.flush();
    }
  }

  cout << endl;
  end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  show_time(begin_ms, end_ms,"TCA calculation");
  spinarrangement_Mathematica_output(randsigma,outfile_spinarr);

  outfile_mlength.close();
  outfile_mcdetails.close();
  outfile_freeenergy.close();
  // outfile_dos.close();
  outfile_spinarr.close();
  return 0;
}